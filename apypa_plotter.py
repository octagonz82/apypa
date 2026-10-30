"""Local, offline interactive viewer for ACAB results. Python standard library only."""
from __future__ import annotations

import argparse
import errno
import io
import json
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import mimetypes
import os
import secrets
from pathlib import Path
import subprocess
import sys
import threading
import time
import urllib.request
import zipfile
from urllib.parse import unquote, urlsplit

from acab_reader import read_fort6, workbook_payload, export_excel


STATIC = Path(__file__).with_name('plotter') / 'static'


def plot_data(result) -> dict:
    """Adapt the reader's authoritative Excel view, without renormalizing values."""
    payload = workbook_payload(result)
    families = {}
    for spec in payload['sheets']:
        matrix = spec['matrix']
        family = spec['family']
        data = {'name': spec['name'], 'unit': spec['unit'],
                'normalization': spec['normalization'], 'samples': [], 'series': []}
        if family == 'gamma':
            data.update(bins=[], values=[], totals=[])
        if matrix and family == 'gamma':
            upper = matrix[0][4:-2]
            boundaries = sorted(result.gamma_boundaries_mev)
            data['bins'] = [{'lower': a, 'upper': b} for a, b in zip(boundaries, boundaries[1:])]
            if upper != [b['upper'] for b in data['bins']]:
                raise ValueError('The gamma boundaries do not match the result columns.')
            data['values'] = []
            data['totals'] = []
            for row in matrix[1:]:
                data['samples'].append({'set': row[0], 'zone': row[1], 'time': row[2], 'label': row[3]})
                data['values'].append(row[4:-2])
                data['totals'].append(row[-2])
        elif matrix:
            data['samples'] = [{'set': s, 'zone': z, 'time': t, 'label': label}
                               for s, z, t, label in zip(matrix[0][1:], matrix[1][1:],
                                                         matrix[2][1:], matrix[3][1:])]
            data['series'] = [{'name': row[0], 'values': row[1:]} for row in matrix[4:]]
        families[family] = data
    return {'cell': result.cell, 'source': f'{result.cell}/fort.6', 'families': families,
            'warnings': result.warnings, 'finished': result.finished,
            'volumes': sorted({t.volume_cm3 for t in result.tables if t.volume_cm3 is not None})}


class ResultStore:
    """Discover cells and cache reads until the source file's fingerprint changes."""
    def __init__(self, root: Path):
        self.root = root.resolve()
        self._lock = threading.Lock()
        self._paths = {}
        self._locks = {}
        self._cache = {}

    def discover(self) -> list[dict]:
        files = ([self.root/'fort.6'] if (self.root/'fort.6').is_file()
                 else sorted(self.root.glob('*/fort.6')))
        paths = {p.parent.name: p for p in files if p.is_file()}
        with self._lock:
            self._paths = paths
            for cell in paths:
                self._locks.setdefault(cell, threading.Lock())
            for cell in set(self._cache) - set(paths):
                del self._cache[cell]
        return [{'id': cell, 'bytes': p.stat().st_size} for cell, p in paths.items()]

    def load(self, cell: str) -> dict:
        with self._lock:
            if cell not in self._paths:
                raise KeyError(f'Cell not available: {cell}')
            path, lock = self._paths[cell], self._locks[cell]
        with lock:
            stat = path.stat()
            fingerprint = (str(path.resolve()), stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino, stat.st_size)
            cached = self._cache.get(cell)
            if cached and cached[0] == fingerprint:
                return cached[1]
            data = plot_data(read_fort6(path))
            after = path.stat()
            if (after.st_mtime_ns, after.st_ctime_ns, after.st_ino, after.st_size) != (
                    stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino, stat.st_size):
                raise ValueError(f'{cell}/fort.6 changed while being read; refresh the data.')
            data['modified_ns'] = stat.st_mtime_ns
            with self._lock:
                if self._paths.get(cell) == path:
                    self._cache[cell] = (fingerprint, data)
            return data

    def change_directory(self, directory):
        candidate = ResultStore(Path(directory))
        cells = candidate.discover()
        if not cells:
            raise ValueError('The directory must contain fort.6 or subdirectories containing fort.6.')
        with self._lock:
            self.root = candidate.root
            self._paths, self._locks, self._cache = candidate._paths, candidate._locks, {}
        return cells


class DesktopActions:
    """Local folder selection and exports, with explicit overwrite selection."""
    def __init__(self, store):
        self.store = store
        self.export_lock = threading.Lock()
        self.downloads = {}

    def choose_directory(self):
        from apypa_launcher import launcher_command
        result = subprocess.run(launcher_command('--pick-directory', str(self.store.root)),
                                capture_output=True, text=True, timeout=300)
        try:
            selected = json.loads(result.stdout)
        except ValueError:
            raise ValueError('Could not open the directory picker. Enter the path in the Data directory field.') from None
        if result.returncode or 'error' in selected:
            raise ValueError('Could not open the directory picker. Enter the path in the Data directory field.')
        return selected

    def export(self, cells, overwrite):
        if not isinstance(cells, list) or not cells or any(not isinstance(c, str) for c in cells):
            raise ValueError('Select at least one cell.')
        if type(overwrite) is not bool:
            raise ValueError('The overwrite option must be a boolean.')
        if not self.export_lock.acquire(blocking=False):
            raise ValueError('An export is already in progress.')
        try:
            with self.store._lock:
                root = self.store.root
                paths = [self.store._paths[c] for c in dict.fromkeys(cells)]
            output = root/'results'
            targets = [output/f'{p.parent.name}.xlsx' for p in paths]
            existing = [p.name for p in targets if p.exists()]
            if existing and not overwrite:
                raise ValueError('Already exist: ' + ', '.join(existing) + '. Enable Overwrite existing to replace them.')
            files, errors = [], []
            for path, target in zip(paths, targets):
                try:
                    export_excel(read_fort6(path), target, overwrite=overwrite)
                    files.append(target)
                except (OSError, ValueError, RuntimeError) as exc:
                    errors.append(f'{path.parent.name}: {exc}')
            response = {'directory': str(output), 'files': [p.name for p in files], 'errors': errors}
            if files:
                archive = io.BytesIO()
                with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as zipped:
                    for file in files:
                        zipped.write(file, file.name)
                token = secrets.token_urlsafe(24)
                self.downloads[token] = archive.getvalue()
                while len(self.downloads) > 3:
                    del self.downloads[next(iter(self.downloads))]
                response['download'] = f'/api/download/{token}'
            return response
        finally:
            self.export_lock.release()


class BrowserSessions:
    """Track live browser connections, allowing a short grace for reloads."""
    def __init__(self, grace=3, timeout=300):
        self.token = secrets.token_urlsafe(24)
        self.grace, self.timeout = grace, timeout
        self.pages = {}
        self.seen = False
        self.idle_since = None
        self.lock = threading.Lock()
        self.closed = threading.Event()

    def update(self, page, sequence, closing, now=None):
        now = time.monotonic() if now is None else now
        with self.lock:
            old = self.pages.get(page)
            if old is None or sequence > old[0]:
                self.pages[page] = (sequence, now, closing)
                self.seen = True
                if not closing:
                    self.idle_since = None

    def should_stop(self, now=None):
        now = time.monotonic() if now is None else now
        with self.lock:
            self.pages = {k: v for k, v in self.pages.items() if now - v[1] < self.timeout}
            if not self.seen or any(not v[2] for v in self.pages.values()):
                self.idle_since = None
                return False
            if self.idle_since is None:
                self.idle_since = now
            return now - self.idle_since >= self.grace


def make_handler(store: ResultStore, sessions=None, actions=None):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt, *args):
            # Only report failed requests; normal polling/assets stay quiet.
            if args and str(args[1] if len(args) > 1 else '').startswith(('4', '5')):
                super().log_message(fmt, *args)

        def send_bytes(self, body, kind, status=HTTPStatus.OK, vendor=False, download=None):
            self.send_response(status)
            self.send_header('Content-Type', kind)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Cache-Control', 'public, max-age=86400' if vendor else 'no-store')
            if download:
                self.send_header('Content-Disposition', f'attachment; filename="{download}"')
            self.end_headers()
            self.wfile.write(body)

        def json(self, data, status=HTTPStatus.OK):
            body = json.dumps(data, ensure_ascii=False, allow_nan=False,
                              separators=(',', ':')).encode('utf-8')
            self.send_bytes(body, 'application/json; charset=utf-8', status)

        def do_GET(self):
            path = unquote(urlsplit(self.path).path)
            try:
                if sessions and path == f'/api/viewer/{sessions.token}':
                    return self.watch_connection()
                if path == '/api/cells':
                    return self.json({'app': 'apypa-plotter', 'root': str(store.root),
                                      'session_token': sessions.token if sessions else None,
                                      'cells': store.discover(), 'pid': os.getpid()})
                if path.startswith('/api/cell/'):
                    cell = path.removeprefix('/api/cell/')
                    return self.json(store.load(cell))
                if actions and path.startswith('/api/download/'):
                    body = actions.downloads[path.removeprefix('/api/download/')]
                    return self.send_bytes(body, 'application/zip', download='APYPA_results.zip')
                if path in ('/', '/index.html'):
                    target = STATIC/'index.html'
                elif path.startswith('/static/'):
                    target = (STATIC/path.removeprefix('/static/')).resolve()
                    if not target.is_relative_to(STATIC.resolve()):
                        raise KeyError('Path not available')
                else:
                    raise KeyError('Path not available')
                if not target.is_file():
                    raise KeyError('File not available')
                kind = mimetypes.guess_type(target)[0] or 'application/octet-stream'
                if kind.startswith('text/') or kind in ('application/javascript',):
                    kind += '; charset=utf-8'
                self.send_bytes(target.read_bytes(), kind, vendor='vendor' in target.parts)
            except KeyError as exc:
                self.json({'error': str(exc).strip("'")}, HTTPStatus.NOT_FOUND)
            except (OSError, ValueError) as exc:
                self.json({'error': str(exc)}, HTTPStatus.UNPROCESSABLE_ENTITY)

        def do_POST(self):
            if not sessions or not actions or self.headers.get('X-APYPA-Token') != sessions.token:
                return self.json({'error': 'Invalid session. Reload the viewer.'}, HTTPStatus.FORBIDDEN)
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if not 0 < length <= 16384:
                    raise ValueError('Invalid request.')
                data = json.loads(self.rfile.read(length))
                if not isinstance(data, dict):
                    raise ValueError('Invalid request.')
                if self.path == '/api/project/pick':
                    return self.json(actions.choose_directory())
                if self.path == '/api/project':
                    store.change_directory(data['directory'])
                    return self.json({'directory': str(store.root)})
                if self.path == '/api/export':
                    return self.json(actions.export(data['cells'], data.get('overwrite', False)))
                return self.json({'error': 'Path not available'}, HTTPStatus.NOT_FOUND)
            except (OSError, ValueError, TypeError, KeyError, subprocess.TimeoutExpired) as exc:
                self.json({'error': str(exc)}, HTTPStatus.UNPROCESSABLE_ENTITY)

        def watch_connection(self):
            page, sequence = secrets.token_urlsafe(24), 0
            self.send_response(HTTPStatus.OK)
            self.send_header('Content-Type', 'text/event-stream')
            self.send_header('Cache-Control', 'no-store')
            self.send_header('Connection', 'close')
            self.end_headers()
            sessions.update(page, sequence, False)
            try:
                while not sessions.closed.is_set():
                    self.wfile.write(b': connected\n\n')
                    self.wfile.flush()
                    sequence += 1
                    sessions.update(page, sequence, False)
                    sessions.closed.wait(1)
            except OSError:
                pass  # Closing a tab breaks its live connection.
            finally:
                sessions.update(page, sequence + 1, True)
                self.close_connection = True

    return Handler


def existing_viewer(port: int, directory: Path):
    """Only reuse an APYPA server serving the requested data directory."""
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/cells', timeout=1) as response:
            data = json.load(response)
        if (data.get('app') == 'apypa-plotter'
                and data.get('root') == str(directory.resolve())):
            return data
    except (OSError, ValueError, AttributeError):
        pass
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description='Local interactive plotter for ACAB results.')
    parser.add_argument('directory', type=Path, nargs='?', default=Path('.'),
                        help='Directory containing cells, or a single cell directory')
    parser.add_argument('--port', type=int, default=8765)
    parser.add_argument('--open', action='store_true', help='Open the browser automatically')
    parser.add_argument('--background', action='store_true', help='Keep the viewer running in the background; save its log in outputs/plotter')
    args = parser.parse_args(argv)
    if not args.directory.is_dir():
        parser.error('The data directory does not exist.')
    if not (1 <= args.port <= 65535):
        parser.error('The port must be between 1 and 65535.')
    store = ResultStore(args.directory)
    cells = store.discover()
    if not (STATIC/'vendor/plotly-cartesian-4.1.2.min.js').is_file():
        parser.error('Missing plotter/static/vendor/plotly-cartesian-4.1.2.min.js.')
    if existing_viewer(args.port, args.directory):
        print(f'The viewer is already running: http://127.0.0.1:{args.port}', flush=True)
        print('Open that address in your browser.', flush=True)
        if args.open:
            from apypa_launcher import open_browser
            open_browser(f'http://127.0.0.1:{args.port}')
        return 0
    if args.background:
        from apypa_launcher import launcher_command
        output = Path(__file__).resolve().parent/'outputs/plotter'
        output.mkdir(parents=True, exist_ok=True)
        log_path = output/'server.log'
        with log_path.open('ab') as log:
            process = subprocess.Popen(launcher_command(str(args.directory.resolve()), '--port', str(args.port), '--no-browser'),
                                       stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                                       start_new_session=True)
        for _ in range(60):
            if process.poll() is not None:
                parser.exit(1, f'Could not start the viewer. See {log_path}.\n')
            try:
                with urllib.request.urlopen(f'http://127.0.0.1:{args.port}/api/cells', timeout=.5) as response:
                    if json.load(response).get('pid') == process.pid:
                        (output/'server.pid').write_text(str(process.pid), encoding='utf-8')
                        print(f'Viewer running: http://127.0.0.1:{args.port} · PID {process.pid}', flush=True)
                        print(f'Log: {log_path}', flush=True)
                        if args.open:
                            from apypa_launcher import open_browser
                            open_browser(f'http://127.0.0.1:{args.port}')
                        return 0
            except (OSError, ValueError):
                pass
            time.sleep(.1)
        process.terminate()
        parser.exit(1, f'The viewer did not respond during startup. See {log_path}.\n')
    sessions = BrowserSessions()
    try:
        server = ThreadingHTTPServer(('127.0.0.1', args.port), make_handler(store, sessions, DesktopActions(store)))
    except OSError as exc:
        # Another invocation may have started after the initial check.
        if existing_viewer(args.port, args.directory):
            print(f'The viewer is already running: http://127.0.0.1:{args.port}', flush=True)
            return 0
        if args.open and exc.errno == errno.EADDRINUSE:
            server = ThreadingHTTPServer(('127.0.0.1', 0), make_handler(store, sessions, DesktopActions(store)))
            args.port = server.server_address[1]
        else:
            parser.exit(1, f'Could not start the viewer: {exc}. Try a different --port.\n')
    print(f'APYPA · {len(cells)} cells · http://127.0.0.1:{args.port}', flush=True)
    print('Open that address in your browser. It stops when the last tab closes; Ctrl+C also stops it.', flush=True)
    stopped = threading.Event()
    def watch_pages():
        while not stopped.wait(.5):
            if sessions.should_stop():
                print('Last tab closed. Viewer stopped.', flush=True)
                server.shutdown()
                return
    threading.Thread(target=watch_pages, daemon=True).start()
    if args.open:
        from apypa_launcher import open_browser
        threading.Thread(target=open_browser, args=(f'http://127.0.0.1:{args.port}',), daemon=True).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\nViewer stopped.')
    finally:
        stopped.set()
        sessions.closed.set()
        server.server_close()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

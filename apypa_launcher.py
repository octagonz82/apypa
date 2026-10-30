"""Desktop entry point, also used by the relocatable Linux executable."""
from pathlib import Path
import json
import os
import sys
import webbrowser


def launcher_command(*arguments):
    if getattr(sys, 'frozen', False):
        return [sys.executable, *arguments]
    return [sys.executable, str(Path(__file__).resolve()), *arguments]


def pick_directory(initial):
    import tkinter as tk
    from tkinter import filedialog
    root = tk.Tk()
    root.withdraw()
    try:
        return filedialog.askdirectory(parent=root, initialdir=str(initial),
                                       title='Choose a directory containing ACAB cells')
    finally:
        root.destroy()


def open_browser(url):
    # Bundled libraries must not be injected into the external system browser.
    previous = os.environ.get('LD_LIBRARY_PATH')
    if getattr(sys, 'frozen', False):
        original = os.environ.get('LD_LIBRARY_PATH_ORIG')
        if original is None:
            os.environ.pop('LD_LIBRARY_PATH', None)
        else:
            os.environ['LD_LIBRARY_PATH'] = original
    try:
        return webbrowser.open(url)
    finally:
        if previous is None:
            os.environ.pop('LD_LIBRARY_PATH', None)
        else:
            os.environ['LD_LIBRARY_PATH'] = previous


def default_directory():
    base = Path(sys.executable).parent if getattr(sys, 'frozen', False) else Path(__file__).resolve().parent
    for candidate in (base.parent, base, Path.cwd()):
        if (candidate/'fort.6').is_file() or any(candidate.glob('*/fort.6')):
            return candidate
    return base.parent


def main(arguments=None):
    arguments = list(sys.argv[1:] if arguments is None else arguments)
    if arguments == ['--diagnostics']:
        import openpyxl
        import tkinter as tk
        root = tk.Tk()
        root.withdraw()
        try:
            print(json.dumps({'python': sys.version.split()[0], 'frozen': bool(getattr(sys, 'frozen', False)),
                              'openpyxl': openpyxl.__version__, 'tk': root.tk.call('info', 'patchlevel')}))
        finally:
            root.destroy()
        return 0
    if arguments[:1] == ['--pick-directory']:
        try:
            print(json.dumps({'directory': pick_directory(arguments[1] if len(arguments) > 1 else Path.home())}))
            return 0
        except Exception as exc:
            print(json.dumps({'error': str(exc)}))
            return 1
    from apypa_plotter import main as serve
    if not arguments or arguments[0].startswith('--'):
        arguments.insert(0, str(default_directory()))
    if '--no-browser' in arguments:
        arguments.remove('--no-browser')
    else:
        arguments.append('--open')
    return serve(arguments)


if __name__ == '__main__':
    raise SystemExit(main())

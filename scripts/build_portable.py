"""Build a relocatable Linux folder using the current virtual environment."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=root/'dist')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='apypa-build-') as directory:
        temporary = Path(directory)
        subprocess.run([sys.executable, '-m', 'PyInstaller', '--noconfirm', '--onedir',
                        '--name', 'APYPA', '--add-data', f'{root / "plotter"}:plotter',
                        '--hidden-import', 'excel_portable', '--collect-all', 'openpyxl',
                        '--distpath', str(output), '--workpath', str(temporary/'work'),
                        '--specpath', str(temporary), str(root/'apypa_launcher.py')],
                       cwd=root, check=True)
    destination = output/'APYPA'
    for name in ('Start_APYPA.sh', 'README_PORTABLE.md', 'THIRD_PARTY_NOTICES.md', 'LICENSE'):
        shutil.copy2(root/name, destination/name)
    # Keep licenses delivered by installed Python distributions with the bundle.
    import importlib.metadata
    licenses = destination/'licenses'
    for package in ('openpyxl', 'et-xmlfile', 'pyinstaller'):
        distribution = importlib.metadata.distribution(package)
        for item in distribution.files or []:
            if any(word in str(item).lower() for word in ('license', 'licence', 'copying')) and str(item).endswith(('.txt', '.rst', '.md', 'COPYING', 'LICENSE')):
                target = licenses/package/Path(str(item)).name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(distribution.locate_file(item), target)
    python_license = Path(sys.base_prefix)/'LICENSE.txt'
    if python_license.is_file():
        licenses.mkdir(exist_ok=True)
        shutil.copyfile(python_license, licenses/'Python-LICENSE.txt')
    collected = {p.name for p in (destination/'_internal').glob('*.so*')}
    # Conda provides exact package provenance for the bundled shared libraries.
    for record in (Path(sys.base_prefix)/'conda-meta').glob('*.json'):
        metadata = json.loads(record.read_text())
        if not any(Path(name).name in collected for name in metadata.get('files', [])):
            continue
        source = Path(metadata.get('link', {}).get('source', ''))/'info/licenses'
        if source.is_dir():
            shutil.copytree(source, licenses/f"{metadata['name']}-{metadata['version']}", dirs_exist_ok=True)
    # Notices for X11 libraries collected from common Debian/Ubuntu hosts.
    for package in ('libx11-6', 'libxau6', 'libxdmcp6', 'libxcb1', 'libbsd0', 'libmd0'):
        notice = Path('/usr/share/doc')/package/'copyright'
        if notice.is_file():
            target = licenses/package/'copyright'
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(notice, target)
    print(f'Portable package created: {destination}')


if __name__ == '__main__':
    main()

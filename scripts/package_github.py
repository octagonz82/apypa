"""Prepare a source-only GitHub handoff using an explicit file allowlist."""
import argparse
from pathlib import Path
import shutil
import zipfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=root/'outputs/APYPA_GitHub')
    args = parser.parse_args()
    destination = args.output.resolve()
    if destination.exists():
        raise FileExistsError(f'The delivery directory already exists: {destination}')
    files = ['acab_reader.py', 'apypa_plotter.py', 'apypa_launcher.py',
             'excel_portable.py', 'excel_export.mjs', 'requirements.txt',
             'requirements-build.txt', '.gitignore', 'Install.sh', 'Start_APYPA.sh',
             'README_PORTABLE.md', 'THIRD_PARTY_NOTICES.md', 'LICENSE']
    for name in files:
        target = destination/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root/name, target)
    shutil.copy2(root/'README.md', destination/'README.md')
    for folder in ('plotter', 'tests', 'scripts', '.github'):
        shutil.copytree(root/folder, destination/folder,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.spec'))
    with zipfile.ZipFile(destination.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(destination.rglob('*')):
            if path.is_file():
                archive.write(path, str(Path(destination.name)/path.relative_to(destination)))
    print(destination)
    print(destination.with_suffix('.zip'))


if __name__ == '__main__':
    main()

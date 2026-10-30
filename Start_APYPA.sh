#!/bin/sh
set -eu
APYPA_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ -x "$APYPA_DIR/APYPA" ]; then
    exec "$APYPA_DIR/APYPA" "$@"
fi
if [ ! -x "$APYPA_DIR/.venv/bin/python" ]; then
    printf 'Run Install.sh first (it installs dependencies in this directory).\n' >&2
    exit 1
fi
exec "$APYPA_DIR/.venv/bin/python" "$APYPA_DIR/apypa_launcher.py" "$@"

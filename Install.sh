#!/bin/sh
set -eu
APYPA_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 -m venv "$APYPA_DIR/.venv"
"$APYPA_DIR/.venv/bin/python" -m pip install -r "$APYPA_DIR/requirements.txt"
printf '\nReady. Open Start_APYPA.sh.\n'

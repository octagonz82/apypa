# APYPA

Acab PYthon PArser by ESS-Bilbao

Local reader and explorer for ACAB `fort.6` results, with interactive charts
and one Excel workbook per cell. Runs without external services or a Codex connection.

## Using the portable version

Extract the Linux package and copy the **entire APYPA directory** into your project,
alongside the cell directories. Open the **APYPA** executable or
**Start_APYPA.sh**. The browser opens and data near the tool are detected.
Use **Browse…**, or enter a path and click **Load**, to select another project.
The process stops when the last tab closes.

```
My_project/
├── APYPA/                  # Complete portable directory
│   ├── APYPA               # Executable
│   └── _internal/          # Bundled dependencies
├── folder1/fort.6
├── folder2/fort.6
└── results/                # Created when Excel workbooks are generated
```

Select the cells and click **Generate Excel**. You can download the workbooks together
as a ZIP file. **Overwrite existing** is disabled by default. See
[README_PORTABLE.md](README_PORTABLE.md) for Linux compatibility and permissions.

## Running from source

Requires Python 3.10 or later. On Linux:

```bash
chmod +x Install.sh "Start_APYPA.sh"
./Install.sh
./"Start_APYPA.sh"
```

Installation creates `.venv` inside the tool directory and downloads openpyxl from
PyPI. Internet access is no longer needed afterwards. You can also use your own environment:

```bash
python -m pip install -r requirements.txt
python apypa_launcher.py "/path/to/my_project"
```

`apypa_launcher.py` opens the browser automatically; `apypa_plotter.py` starts
the server and displays the link. Both accept `--port 8766`. The launcher
chooses a free port if the usual one is occupied by another project. You can
start with a directory containing no data and select the data from the viewer.
To use just the reader and exporter:

```bash
python acab_reader.py "/path/to/project" --output "/path/to/project/results"
python acab_reader.py "/path/to/project" --output results --overwrite
```

The reader uses the standard library. Excel workbooks use openpyxl, which is
included in the portable package; they do not require Node or the Codex runtime.
`excel_export.mjs` is retained solely for development previews using
artifact-tool; it is not used by the portable application.

## Results and physical conventions

Each Excel workbook contains five worksheets:

| Family | Unit | Rows | Columns |
|---|---|---|---|
| Composition | mol/cm³ | Isotopes | Times |
| Activity | Bq/cm³ | Isotopes | Times |
| Thermal power | W/cm³ | Isotopes | Times |
| Surface gamma dose | Sv/h | Isotopes | Times |
| Gamma spectrum | photons/(cm³·s) | Times | Upper energy boundaries [MeV] |

Composition, activity and thermal power are divided by the printed volume of
each zone. One gram-atom equals one mole. Dose and spectrum retain their
original values; the spectrum is already normalized by volume.

Numeric times are retained in ascending order; INITIAL, RESTART and
SHUTDOWN are excluded from the views. An additional row or column shows the
equivalent duration, using 365-day years and 30-day months. Zero components are
omitted. Times and values that ACAB did not print are not reconstructed.

Isotopes whose visible values are all zero and the SUBTOT row are omitted.
TOTAL is retained as a summary and must not be added to the isotope values.
If there is only one zone, it is omitted from legends and Excel headers.
The original information remains available in `read_fort6()`.

The viewer allows you to compare cells, select the top N isotopes plus TOTAL,
search for isotopes, choose times or gamma groups, adjust colors, increase chart
and legend font sizes, set any X/Y endpoint, use logarithmic scales and export
CSV, PNG or SVG files. Zero values are retained in the data; they are hidden only
when a logarithmic scale prevents plotting them. Maximum: 80 simultaneous curves.
The time filter reduces the observations; X/Y limits only change the chart window.

The server listens exclusively on `127.0.0.1`. The Browse button opens a native
directory picker; if unavailable, you can enter the path. Directory changes and
exports require the local session token. User data and the ACAB manual are not
included in this repository.

## Building the portable package

Build on the target operating system and architecture. On Linux:

```bash
python -m venv .venv-build
.venv-build/bin/python -m pip install -r requirements-build.txt
.venv-build/bin/python scripts/build_portable.py
```

The `dist/APYPA` directory contains the executable and its dependencies. Copy the
entire directory, including the dependencies. The delivered local build was tested on
Ubuntu 24.04 x86_64. The GitHub Actions workflow also generates a Linux artifact;
a separate version must be built and tested for Windows or macOS.

## Checks

```bash
python -m unittest discover -s tests -v
node tests/test_plotter_model.mjs
```

Tests using synthetic data do not require external files. Optional tests using
the ten original cells or previously saved Excel workbooks are skipped when
those data are unavailable. The reader reports truncated tables, invalid numbers,
unknown labels and gamma group mismatches; it does not silently turn these errors
into zeros.

## Components

- `acab_reader.py`: reading and physical transformations.
- `excel_portable.py`: XLSX formatting and export.
- `apypa_plotter.py`: local server and project operations.
- `apypa_launcher.py`: launcher and directory picker.
- `plotter/static/`: interface and local Plotly assets.

CONTRIBUTING

Contributions are highly welcome! If you would like to contribute, it is as easy as forking the repository on GitHub, making your changes, and issuing a pull request. If you have any questions about this process don't hesitate to ask the author.

## License

APYPA is distributed under the GNU General Public License, version 3 (GPLv3).
See [LICENSE](LICENSE) for the complete license text.
When distributing the executable, also provide the corresponding source code
and build scripts in accordance with GPLv3. For downloads, state next to the
executable where that source can be obtained at no additional charge.
Third-party components retain their own licenses and notices, documented in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

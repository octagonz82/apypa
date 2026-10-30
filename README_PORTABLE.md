# Portable APYPA for Linux

1. Copy the entire APYPA directory into your project (alongside the cell directories).
2. Open the `APYPA` executable or run `Start_APYPA.sh`. The browser opens.
3. If no data are detected, use **Browse…**, or enter the path and click **Load**.
4. Select cells and click **Generate Excel**. Workbooks are saved in `results/`
   inside the data directory. They can also be downloaded together as a ZIP file.
5. When the last tab closes, Python exits and releases the port.

You can move this directory. Keep `APYPA`, `_internal` and all other files together.
The executable includes Python, openpyxl, Plotly and the directory picker;
it does not require Codex, Node, a Python installation or an Internet connection.

This package targets Linux x86_64 and was built and tested on Ubuntu 24.04
(glibc 2.39). It is not a universal executable for Windows, macOS or every Linux
distribution. For another system, build the source code on that system.
A graphical session and a browser on the computer are required. If the file manager
does not run programs on double-click, use a terminal: `./APYPA`. If extraction
removes execution permissions, run: `chmod +x APYPA "Start_APYPA.sh"`.
To check the bundled Python, openpyxl and graphical picker, run: `./APYPA --diagnostics`.

Chart limits and preferences are retained during the session. Data are processed
locally on your computer; they are not sent to external services.

The data directory can contain `fort.6` directly or subdirectories named
`cell_name/fort.6`. Each subdirectory name is used for its Excel workbook.
**Overwrite existing** is disabled by default.

The export includes all five families, using the units and normalizations
explained in the source README. INITIAL, RESTART, SHUTDOWN, SUBTOT and isotopes
whose values are always zero are excluded. A single zone is omitted from headers
and legends; it remains available in the internal reader.

## License

APYPA is distributed under the GNU General Public License, version 3 (GPLv3).
See [LICENSE](LICENSE) for the complete license text.
When distributing the executable, also provide the corresponding source code
and build scripts in accordance with GPLv3. For downloads, state next to the
executable where that source can be obtained at no additional charge.
Third-party components retain their own licenses and notices, documented in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

# Third-party components

- Plotly.js Cartesian 4.1.2 and its adapted English locale module: MIT license.
  The original license text is included in `plotter/static/vendor/LICENSE.plotly.txt`
  (in the executable: `_internal/plotter/static/vendor/`).
  The locale module was adapted from the bundled Spanish translation to use
  English labels and formatting.
  https://github.com/plotly/plotly.js
- openpyxl 3.1.5 and et-xmlfile 2.0.0: MIT license.
  https://openpyxl.readthedocs.io/ and https://pypi.org/project/et-xmlfile/
- Python and its standard library: Python Software Foundation license and
  additional component notices. https://docs.python.org/3/license.html
- The directory picker uses Tcl/Tk through tkinter. Tcl/Tk uses permissive
  BSD licenses. https://www.tcl-lang.org/software/tcltk/license.html
- The executable is built with PyInstaller 6.22.3, licensed under the GPL with a
  bootloader exception allowing bundled programs to be distributed under their own license.
  https://pyinstaller.org/en/stable/license.html

The portable package's `licenses/` directory preserves the license texts available
in the Python distributions used during the build. Dynamic libraries are bundled
by PyInstaller; their accompanying component notices are part of the distribution.
APYPA's own code is distributed under GPLv3; see [LICENSE](LICENSE).
The licenses above apply to third-party components and are retained.

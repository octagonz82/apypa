"""Lossless reader for the five main ACAB fort.6 result families.

Extraction uses only the Python standard library. Portable Excel export uses
openpyxl. The optional development preview uses the artifact-tool Node runtime.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, field
from decimal import Decimal
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


FAMILIES = {
    "composition": ("Composition", "Gram-atoms", "NUCLIDE CONCENTRATIONS,  AT.GR."),
    "activity": ("Activity", "Bq", "NUCLIDE RADIOACTIVITY, DISINTEGRATIONS/SEC"),
    "heat": ("Thermal power", "W", "NUCLIDE THERMAL POWER, WATTS"),
    "dose": ("Surface gamma dose", "Sv/h", "SURFACE GAMMA DOSE RATES DUE TO"),
    "gamma": ("Gamma spectrum", "photons/cm³/s", "GROUP PHOTON RELEASE RATES, PHOTONS/CCM/SEC"),
}
TIME_SET = re.compile(r"(\d+)\.\s*TIME SET\s*$")
ZONE = re.compile(r"^\s*ZONE\s+(\d+)\s*$")
VOLUME = re.compile(r"VOLUME OF ZONE\s*[:%]\s*(\S+)")
NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:(?:[EeDd][+-]?\d+)|(?:[+-]\d{2,3}))?")
NUCLIDE = re.compile(r"([A-Za-z]{1,2})\s*(\d+)([A-Za-z]\d*)?")
STATES = {"INITIAL", "RESTART", "SHUTDOWN"}


class ParseError(ValueError):
    """A result cannot be interpreted without losing or inventing data."""


def number(token: str) -> float:
    """Accept Fortran E/D exponents and the E-less exponent form, never default to 0."""
    token = token.strip()
    if not NUMBER.fullmatch(token):
        raise ParseError(f"Invalid number: {token!r}")
    value = token.replace("D", "E").replace("d", "e")
    if "e" not in value.lower():
        value = re.sub(r"(?<=\d)([+-]\d{2,3})$", r"E\1", value)
    result = float(value)
    if not math.isfinite(result):
        raise ParseError(f"Non-finite number: {token!r}")
    if result == 0 and any(c in "123456789" for c in re.split(r"[Ee]", value)[0]):
        raise ParseError(f"Value too small for float64: {token!r}")
    return result


def numbers(text: str, count: int) -> list[float]:
    matches = list(NUMBER.finditer(text))
    remainder = NUMBER.sub("", text)
    if remainder.strip() or len(matches) != count:
        raise ParseError(f"Expected {count} values: {text.strip()!r}")
    return [number(m.group()) for m in matches]


@dataclass(frozen=True)
class TimePoint:
    label: str
    seconds: float | None


@dataclass
class ResultTable:
    family: str
    time_set: int
    zone: int
    volume_cm3: float | None
    line: int
    times: list[TimePoint]
    rows: dict[str, list[float]]
    energies_mev: list[float] = field(default_factory=list)
    dose_description: str = ""


@dataclass
class AcabResult:
    source: Path
    cell: str
    tables: list[ResultTable]
    gamma_boundaries_mev: list[float]
    warnings: list[str]
    finished: bool

    def family(self, name: str) -> list[ResultTable]:
        return [t for t in self.tables if t.family == name]


def _family(line: str) -> str | None:
    normalized = " ".join(line.split()).lstrip("0 ")
    for key, (_, _, heading) in FAMILIES.items():
        if " ".join(heading.split()) in normalized:
            return key
    return None


def _time_header(line: str, family: str) -> list[TimePoint] | None:
    tokens = line.split()
    if not tokens:
        return None
    if family == "gamma":
        if tokens[0] != "(MEV)":
            return None
        tokens = tokens[1:]
    elif tokens[0] == "ISOTOPE":
        tokens = tokens[1:]
    elif tokens[0] not in STATES and not re.fullmatch(r"\S+[Ss]", tokens[0]):
        return None
    if not tokens:
        return None
    times = []
    for token in tokens:
        if token in STATES:
            times.append(TimePoint(token, None))
        elif token.endswith(("S", "s")):
            times.append(TimePoint(token, number(token[:-1])))
        else:
            raise ParseError(f"Unrecognized time: {token!r}")
    return times


def read_fort6(path: str | Path) -> AcabResult:
    """Read all printed samples, including repeated states and zero nuclides.

    Time values are exactly those printed in the result headers, in seconds;
    INITIAL/RESTART/SHUTDOWN remain distinct labels, without invented times.
    """
    path = Path(path).resolve()
    raw = path.read_bytes()
    warnings = []
    if b"\0" in raw:
        warnings.append(f"{raw.count(bytes([0]))} NUL characters replaced with spaces for reading.")
    lines = raw.decode("utf-8").replace("\0", " ").splitlines()
    tables = []
    time_set = None
    zone = None
    volume = None
    nogg = None
    boundaries = []
    i = 0
    while i < len(lines):
        line = lines[i]
        match = TIME_SET.search(line)
        if match:
            time_set = int(match[1])
            volume = None
        match = ZONE.match(line)
        if match:
            zone = int(match[1])
        match = VOLUME.search(line)
        if match:
            volume = number(match[1])
        if "NOGG   NUMBER OF ACTIVATION GAMMA GROUPS" in line:
            nogg = int(line.split()[-1])
        if re.match(r"^\s*BLOCK #2,\s*CARD #6\s*$", line):
            j = i + 1
            while j < len(lines) and lines[j].strip():
                boundaries.extend(number(t) for t in lines[j].split())
                j += 1
        family = _family(line)
        if family is None:
            i += 1
            continue
        first = i
        times = None
        description = ""
        i += 1
        while i < len(lines):
            if TIME_SET.search(lines[i]) or _family(lines[i]):
                break
            match = ZONE.match(lines[i])
            if match:
                zone = int(match[1])
            if family == "dose" and "SV/H" in lines[i]:
                description = lines[i].strip()
            times = _time_header(lines[i], family)
            if times:
                break
            i += 1
        if not times or time_set is None or zone is None:
            raise ParseError(f"{path}:{first+1}: missing times, time set or zone.")
        rows = {}
        energies = []
        i += 1
        while i < len(lines):
            row = lines[i]
            if not row.strip():
                i += 1
                continue
            try:
                if family == "gamma":
                    summary = re.match(r"^\s*0?\s*(TOTAL|MEV/SEC)\s+(.+)$", row)
                    if summary:
                        label, values = summary[1], numbers(summary[2], len(times))
                    else:
                        fields = row.split(maxsplit=1)
                        if len(fields) != 2 or not NUMBER.fullmatch(fields[0]):
                            break
                        energy = number(fields[0])
                        label, values = fields[0], numbers(fields[1], len(times))
                        energies.append(energy)
                else:
                    # ACAB uses a seven-character nuclide field, including M/S states.
                    label = row[:7].replace(" ", "")
                    if label not in {"TOTAL", "SUBTOT"} and not NUCLIDE.fullmatch(label):
                        if not re.fullmatch(r"\d+[A-Za-z]?", label):
                            raise ParseError(f"Unrecognized isotope row: {row!r}")
                        warning = f"Label without a chemical symbol in the source: {label}; preserved without assigning an element."
                        if warning not in warnings:
                            warnings.append(warning)
                    if label not in {"TOTAL", "SUBTOT"}:
                        label = label.capitalize()
                    values = numbers(row[7:], len(times))
                if label in rows:
                    raise ParseError(f"Duplicate row: {label}")
                rows[label] = values
            except ParseError as exc:
                raise ParseError(f"{path}:{i+1}: {exc}") from exc
            i += 1
            if family != "gamma" and label == "TOTAL":
                break
            if family == "gamma" and label == "MEV/SEC":
                break
        if "TOTAL" not in rows:
            raise ParseError(f"{path}:{first+1}: incomplete table, missing TOTAL.")
        if family == "gamma":
            if nogg is not None and len(energies) != nogg:
                raise ParseError(f"{path}:{first+1}: {len(energies)} gamma groups; NOGG={nogg}.")
            if len(set(energies)) != len(energies):
                raise ParseError(f"{path}:{first+1}: duplicate mean energies.")
        tables.append(ResultTable(family, time_set, zone, volume, first+1, times,
                                  rows, energies, description))
    if not tables:
        raise ParseError(f"{path}: no ACAB result tables found.")
    if boundaries and nogg is not None and len(boundaries) != nogg + 1:
        raise ParseError(f"{path}: expected NOGG+1 gamma boundaries.")
    if boundaries and any(a <= b for a, b in zip(boundaries, boundaries[1:])):
        raise ParseError(f"{path}: gamma boundaries are not strictly decreasing.")
    missing = set(FAMILIES) - {t.family for t in tables}
    if missing:
        warnings.append("Families not printed: " + ", ".join(sorted(missing)))
    finished = any("***** JOB FINISHED *****" in l for l in lines)
    if not finished:
        warnings.append("JOB FINISHED not found; the calculation may be incomplete.")
    if any("WASTE DISPOSAL RATINGS" in l for l in lines):
        warnings.append("The file also contains WASTE DISPOSAL RATINGS, outside the five requested worksheets.")
    return AcabResult(path, path.parent.name, tables, boundaries, warnings, finished)


def time_text(seconds: float) -> str:
    """English duration: 365-day years, 30-day months; omit zero components.

    Minutes are intentionally not used: the remainder after hours is in seconds.
    The numeric time remains the authoritative value in the adjacent row/column.
    """
    if not math.isfinite(seconds) or seconds < 0:
        raise ValueError("Time must be finite and non-negative.")
    remaining = Decimal(str(seconds))
    parts = []
    units = [(31536000, 'year', 'years'), (2592000, 'month', 'months'),
             (86400, 'day', 'days'), (3600, 'hour', 'hours')]
    for scale, singular, plural in units:
        count = int(remaining // scale)
        remaining -= count * scale
        if count:
            parts.append(f"{count} {singular if count == 1 else plural}")
    if remaining:
        text = format(remaining, 'f')
        if '.' in text:
            text = text.rstrip('0').rstrip('.')
        parts.append(f"{text} {'second' if remaining == 1 else 'seconds'}")
    return ', '.join(parts[:-1]) + ' and ' + parts[-1] if len(parts) > 1 else (parts[0] if parts else '')


def workbook_payload(result: AcabResult) -> dict:
    """Build physically labelled matrices; omit all-zero nuclides in Excel views.

    Absent observations are null/blank, never fabricated zeros. Excel shows only
    numeric times, sorted ascending; the reader still preserves all source states.
    Only extensive composition, activity and heat are divided by volume. Dose and
    the already-volumetric gamma source retain their original values.
    Gamma columns use ascending upper group boundaries from the input.
    """
    sheets = []
    normalized_units = {"composition": "mol/cm³", "activity": "Bq/cm³",
                        "heat": "W/cm³", "dose": "Sv/h",
                        "gamma": "photons/(cm³·s)"}
    for family, (name, unit, _) in FAMILIES.items():
        tables = result.family(family)
        divide_volume = family in {"composition", "activity", "heat"}
        sheet = {"name": name, "family": family, "unit": normalized_units[family],
                 "source_unit": unit,
                 "normalization": ("printed value / zone volume [cm³]" if divide_volume
                                   else "original values, without additional volume division"),
                 "source": str(result.source), "cell": result.cell}
        if not tables:
            sheet.update(matrix=[], retained=0, omitted=0)
            sheets.append(sheet)
            continue
        samples = sorted(
            [(t, j, p.seconds) for t in tables for j, p in enumerate(t.times)
             if p.seconds is not None], key=lambda sample: sample[2])
        if not samples:
            sheet.update(matrix=[], retained=0, omitted=0)
            sheets.append(sheet)
            continue
        for table in tables if divide_volume else []:
            if table.volume_cm3 is None or not math.isfinite(table.volume_cm3) or table.volume_cm3 <= 0:
                raise ParseError(f"{result.source}:{table.line}: a positive volume is required for normalization.")
        if family == "gamma":
            boundaries = result.gamma_boundaries_mev
            if len(boundaries) < 2 or any(not math.isfinite(e) or e < 0 for e in boundaries) or any(
                    a <= b for a, b in zip(boundaries, boundaries[1:])):
                raise ParseError(f"{result.source}: valid gamma boundaries are required to display bin upper edges.")
            energies = sorted(boundaries[:-1])
            for table in tables:
                if len(table.energies_mev) != len(energies) or any(
                        not math.isclose(mean, (upper+lower)/2, rel_tol=0.005, abs_tol=1e-12)
                        for mean, upper, lower in zip(table.energies_mev, boundaries, boundaries[1:])):
                    raise ParseError(f"{result.source}:{table.line}: the printed groups do not match the gamma boundaries.")
            matrix = [["Time set", "Zone", "Time [s]", "Equivalent duration", *energies,
                       "TOTAL [photons/(cm³·s)]", "Gamma power [MeV/(cm³·s)]"]]
            for table, j, seconds in samples:
                by_energy = {e: v for e, v in zip(boundaries[:-1],
                             (v for k, v in table.rows.items() if k not in {"TOTAL", "MEV/SEC"}))}
                matrix.append([table.time_set, table.zone, seconds, time_text(seconds),
                                   *[by_energy[e][j] for e in energies],
                                   table.rows["TOTAL"][j],
                                   table.rows["MEV/SEC"][j] if "MEV/SEC" in table.rows else None])
            sheet["energy_header"] = "Bin upper edge [MeV], ascending order"
            sheet.update(matrix=matrix, retained=len(energies), omitted=0)
        else:
            nuclides = list(dict.fromkeys(k for t in tables for k in t.rows if k not in {"TOTAL", "SUBTOT"}))
            keep = [k for k in nuclides if any(k in t.rows and t.rows[k][j] != 0
                                             for t, j, _ in samples)]
            # SUBTOT remains in the raw reader for checks, outside public views.
            # Retain TOTAL as printed, even when zero; it is not an isotope.
            labels = keep + (["TOTAL"] if any("TOTAL" in t.rows for t in tables) else [])
            matrix = [
                ["Time set", *[t.time_set for t, _, _ in samples]],
                ["Zone", *[t.zone for t, _, _ in samples]],
                ["Isotope / time [s]", *[seconds for _, _, seconds in samples]],
                ["Equivalent duration", *[time_text(seconds) for _, _, seconds in samples]],
            ]
            matrix.extend([k, *[t.rows[k][j] / (t.volume_cm3 if divide_volume else 1) if k in t.rows else None
                                for t, j, _ in samples]] for k in labels)
            sheet.update(matrix=matrix, retained=len(keep), omitted=len(nuclides)-len(keep))
        sheet["volumes"] = list(dict.fromkeys(t.volume_cm3 for t in tables))
        sheet["dose_descriptions"] = list(dict.fromkeys(t.dose_description for t in tables if t.dose_description))
        sheet["table_sources"] = [{"line": t.line, "time_set": t.time_set, "zone": t.zone,
                                   "volume_cm3": t.volume_cm3, "times": [p.label for p in t.times]}
                                  for t in tables]
        sheets.append(sheet)
    return {"cell": result.cell, "sheets": sheets,
            "zones": sorted({t.zone for t in result.tables}),
            "gamma_boundaries_mev": result.gamma_boundaries_mev,
            "warnings": result.warnings, "finished": result.finished}


def _runtime() -> tuple[str, Path]:
    bundle = Path.home() / ".cache/codex-runtimes/codex-primary-runtime/dependencies/node"
    node = os.environ.get("APYPA_NODE") or str(bundle / "bin/node")
    modules = Path(os.environ.get("APYPA_NODE_MODULES", str(bundle / "node_modules")))
    if not Path(node).is_file():
        node = shutil.which("node") or ""
    if not node or not (modules / "@oai/artifact-tool").is_dir():
        raise RuntimeError("The @oai/artifact-tool Excel engine is missing. Set APYPA_NODE and APYPA_NODE_MODULES. Reading with read_fort6() does not require Node.")
    return node, modules


def export_excel(result: AcabResult, destination: str | Path, *, overwrite: bool = False,
                 preview_dir: str | Path | None = None) -> Path:
    """Export five worksheets, with transactional replacement of the workbook."""
    if preview_dir is None:
        from excel_portable import write_workbook
        return write_workbook(workbook_payload(result), destination, overwrite=overwrite)
    destination = Path(destination).resolve()
    if destination.exists() and not overwrite:
        raise FileExistsError(f"{destination} already exists; use --overwrite to replace it.")
    node, modules = _runtime()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="apypa-") as temp:
        temp = Path(temp)
        (temp / "node_modules").symlink_to(modules, target_is_directory=True)
        shutil.copyfile(Path(__file__).with_name("excel_export.mjs"), temp / "excel_export.mjs")
        payload = workbook_payload(result)
        payload["output"] = str(temp / "result.xlsx")
        if preview_dir is not None:
            payload["preview_dir"] = str(Path(preview_dir).resolve())
        (temp / "input.json").write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False), encoding="utf-8")
        subprocess.run([node, str(temp / "excel_export.mjs"), str(temp / "input.json")], check=True)
        # Stage on the destination filesystem to make replacement atomic.
        with tempfile.NamedTemporaryFile(dir=destination.parent, suffix=".xlsx", delete=False) as staged:
            stage_path = Path(staged.name)
        try:
            shutil.copyfile(temp / "result.xlsx", stage_path)
            if overwrite:
                os.replace(stage_path, destination)
            else:
                os.link(stage_path, destination)  # Do not overwrite a concurrent export.
        finally:
            stage_path.unlink(missing_ok=True)
    return destination


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read ACAB fort.6 files and create a five-sheet Excel workbook per cell.")
    parser.add_argument("input", nargs="?", default=".", type=Path, help="fort.6, a cell directory or a root directory")
    parser.add_argument("--output", type=Path, default=Path("outputs/excel"), help="Output directory")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing Excel workbooks")
    parser.add_argument("--preview-dir", type=Path, help="Save PNG previews of each worksheet for review")
    args = parser.parse_args(argv)
    if args.input.is_file():
        paths = [args.input]
    elif (args.input / "fort.6").is_file():
        paths = [args.input / "fort.6"]
    else:
        paths = sorted(args.input.glob("*/fort.6"))
    if not paths:
        parser.error(f"No fort.6 files in {args.input}")
    failures = 0
    for path in paths:
        try:
            result = read_fort6(path)
            target = args.output / f"{result.cell}.xlsx"
            preview = args.preview_dir / result.cell if args.preview_dir else None
            export_excel(result, target, overwrite=args.overwrite, preview_dir=preview)
            report = workbook_payload(result)
            (args.output / f"{result.cell}.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
            counts = Counter(t.family for t in result.tables)
            print(f"{result.cell}: {target} ({dict(counts)})", flush=True)
            for warning in result.warnings:
                print(f"  Warning: {warning}", flush=True)
        except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
            failures += 1
            print(f"ERROR {path}: {exc}", flush=True)
    return int(failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())

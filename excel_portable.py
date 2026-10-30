"""Portable XLSX rendering; the scientific transformations live in acab_reader."""
from pathlib import Path
import os
import tempfile


def visible_matrix(payload, spec):
    matrix = spec['matrix']
    if len(payload['zones']) > 1:
        return matrix
    if spec['family'] == 'gamma':
        return [[value for index, value in enumerate(row) if index != 1] for row in matrix]
    return [row for index, row in enumerate(matrix) if index != 1]


def write_workbook(payload, destination, *, overwrite=False):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    destination = Path(destination).resolve()
    if destination.exists() and not overwrite:
        raise FileExistsError(f'{destination} already exists; enable overwrite to replace it.')
    workbook = Workbook()
    workbook.remove(workbook.active)
    show_zone = len(payload['zones']) > 1
    body_font = Font(name='Arial', size=10)
    header_font = Font(name='Arial', size=10, bold=True, color='FFFFFF')
    total_font = Font(name='Arial', size=10, bold=True)
    header_fill = PatternFill('solid', fgColor='233B55')
    for spec in payload['sheets']:
        sheet = workbook.create_sheet(spec['name'])
        sheet.sheet_view.showGridLines = False
        gamma = spec['family'] == 'gamma'
        divided = spec['family'] in ('composition', 'activity', 'heat')
        sheet['A1'] = f"Cell {payload['cell']} · {spec['name']}"
        sheet['A1'].font = Font(name='Arial', size=14, bold=True)
        sheet['A2'] = f"Source: {spec['source']}"
        sheet['A3'] = f"Units: {spec['unit']}" + (f" · {spec.get('energy_header', 'Bin upper edge [MeV]')}" if gamma else '')
        sheet['A4'] = (f"{' / '.join(spec.get('dose_descriptions', []))} · Original dose in Sv/h."
                       if spec['family'] == 'dose' else
                       'Volumetric source per group; original ACAB values.' if gamma else
                       f"Isotopes retained: {spec['retained']} · Omitted as all zero: {spec['omitted']}")
        volume = ', '.join(str(v) for v in spec.get('volumes', []) if v is not None)
        sheet['A5'] = (f"Volume [cm³]: {volume or 'not printed'} · "
                       f"{'Results / volume' if divided else 'No additional volume division'}"
                       ' · Year: 365 days; month: 30 days.')
        for row in range(2, 6):
            sheet.cell(row, 1).font = Font(name='Arial', size=10, color='526174')
        matrix = visible_matrix(payload, spec)
        if not matrix:
            sheet['A7'] = 'No results with numeric times in this family.'
            continue
        columns = len(matrix[0])
        if len(matrix) + 5 > 1048576 or columns > 16384:
            raise ValueError(f"Worksheet {spec['name']} exceeds Excel limits.")
        if any(len(row) != columns for row in matrix):
            raise ValueError(f"Irregular matrix: {spec['name']}")
        headers = 1 if gamma else 4 if show_zone else 3
        for c in range(1, columns + 1):
            sheet.column_dimensions[get_column_letter(c)].width = 16
        sheet.column_dimensions['A'].width = 20 if gamma else 28
        for i, row in enumerate(matrix):
            number = i + 6
            sheet.row_dimensions[number].height = 34 if gamma and i == 0 else 24 if i < headers else 32 if gamma else 19
            for j, value in enumerate(row, 1):
                cell = sheet.cell(number, j, value)
                if isinstance(value, str):
                    cell.data_type = 's'
                cell.font = body_font
                cell.alignment = Alignment(vertical='center')
                if i < headers:
                    cell.fill, cell.font = header_fill, header_font
                    cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
                    cell.border = Border(right=Side(style='thin', color='FFFFFF'))
                if not gamma and i >= headers and row[0] == 'TOTAL':
                    cell.fill = PatternFill('solid', fgColor='E5ECF3')
                    cell.font = total_font
                if isinstance(value, (int, float)) and (
                        gamma and (j == (3 if show_zone else 2) and i > 0 or j >= (5 if show_zone else 4))
                        or not gamma and j > 1 and (i == headers - 2 or i >= headers)):
                    cell.number_format = '0.000E+00'
        if gamma:
            time_col = 3 if show_zone else 2
            if show_zone:
                sheet.column_dimensions['B'].width = 9
            sheet.column_dimensions[get_column_letter(time_col)].width = 26
            sheet.column_dimensions[get_column_letter(time_col + 1)].width = 55
            for row in sheet.iter_rows(min_row=7, min_col=time_col + 1, max_col=time_col + 1):
                row[0].alignment = Alignment(vertical='center', wrap_text=True)
            for c in (columns - 1, columns):
                sheet.column_dimensions[get_column_letter(c)].width = 32
            sheet.freeze_panes = 'D7' if show_zone else 'C7'
        else:
            sheet.row_dimensions[5 + headers].height = 72
            sheet.freeze_panes = 'B10' if show_zone else 'B9'
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, suffix='.xlsx', delete=False) as staged:
        stage = Path(staged.name)
    try:
        workbook.save(stage)
        if overwrite:
            os.replace(stage, destination)
        else:
            os.link(stage, destination)
    finally:
        stage.unlink(missing_ok=True)
        workbook.close()
    return destination

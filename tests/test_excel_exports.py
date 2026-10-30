"""Optional end-to-end checks of saved Excel files (openpyxl is read-only here)."""
import json
import math
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET
from zipfile import ZipFile
from acab_reader import time_text

try:
    import openpyxl
except ImportError:
    openpyxl = None

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'outputs/excel'
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


@unittest.skipUnless(openpyxl is not None, 'openpyxl is only required to verify saved Excel workbooks')
class SavedExcelTests(unittest.TestCase):
    def test_every_exported_value_and_layout(self):
        paths = sorted(OUTPUT.glob('*.xlsx'))
        if not paths:
            self.skipTest('First generate Excel workbooks with python3 acab_reader.py')
        self.assertEqual([p.stem for p in paths], [str(n) for n in range(6000, 6010)])
        for path in paths:
            with self.subTest(cell=path.stem):
                spec = json.loads(path.with_suffix('.json').read_text())
                workbook = openpyxl.load_workbook(path, read_only=True, data_only=False)
                self.assertEqual(workbook.sheetnames, [s['name'] for s in spec['sheets']])
                with ZipFile(path) as archive:
                    for i, (sheet, expected) in enumerate(zip(workbook, spec['sheets']), 1):
                        source_matrix = expected['matrix']
                        show_zone = len(spec['zones']) > 1
                        matrix = (source_matrix if show_zone else
                                  [[v for j, v in enumerate(row) if j != 1] for row in source_matrix]
                                  if expected['family'] == 'gamma' else
                                  [row for j, row in enumerate(source_matrix) if j != 1])
                        if not show_zone:
                            self.assertNotIn('Zone', matrix[0] if expected['family'] == 'gamma'
                                             else [row[0] for row in matrix])
                        times = ([row[2] for row in source_matrix[1:]] if expected['family'] == 'gamma'
                                 else source_matrix[2][1:])
                        self.assertTrue(all(isinstance(t, (int, float)) for t in times))
                        self.assertEqual(times, sorted(times))
                        self.assertEqual(len(times), 60)
                        if expected['family'] == 'gamma':
                            self.assertEqual(source_matrix[0][4:-2], sorted(spec['gamma_boundaries_mev'][:-1]))
                            self.assertEqual([row[3] for row in source_matrix[1:]], [time_text(t) for t in times])
                        else:
                            self.assertNotIn('SUBTOT', [row[0] for row in source_matrix[4:]])
                            self.assertEqual(source_matrix[3][1:], [time_text(t) for t in times])
                        # OOXML dimensions are optional; streaming readers must calculate them.
                        sheet.calculate_dimension(force=True)
                        self.assertEqual(sheet.max_row, len(matrix)+5)
                        self.assertEqual(sheet.max_column, len(matrix[0]))
                        self.assertEqual(sheet['A3'].value.split(' · ')[0], f"Units: {expected['unit']}")
                        actual = sheet.iter_rows(min_row=6, max_row=sheet.max_row,
                                                 max_col=sheet.max_column)
                        for row_index, (row, expected_row) in enumerate(zip(actual, matrix), 6):
                            for cell, value in zip(row, expected_row):
                                if isinstance(value, (int, float)):
                                    self.assertIsInstance(cell.value, (int, float))
                                    self.assertTrue(math.isclose(cell.value, value, rel_tol=3e-15, abs_tol=0),
                                                    f'{path.stem}/{sheet.title}/{cell.coordinate}')
                                    self.assertEqual(cell.data_type, 'n')
                                else:
                                    self.assertEqual(cell.value, value,
                                                     f'{path.stem}/{sheet.title}/{cell.coordinate}')
                                self.assertNotEqual(cell.data_type, 'e')
                        xml = ET.fromstring(archive.read(f'xl/worksheets/sheet{i}.xml'))
                        pane = xml.find('.//s:pane', NS)
                        expected_pane = ('D7' if show_zone else 'C7') if expected['family'] == 'gamma' else ('B10' if show_zone else 'B9')
                        self.assertIsNotNone(pane)
                        self.assertEqual(pane.attrib['topLeftCell'], expected_pane)
                workbook.close()


if __name__ == '__main__':
    unittest.main()

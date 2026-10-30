import json
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile

import openpyxl

from acab_reader import read_fort6, workbook_payload, export_excel
from apypa_plotter import ResultStore, DesktopActions, BrowserSessions
from excel_portable import visible_matrix

BODY = ('1Test 6. TIME SET\n VOLUME OF ZONE : 10 CCM\n'
        '0 NUCLIDE RADIOACTIVITY, DISINTEGRATIONS/SEC\n ZONE 1\n'
        ' INITIAL RESTART 1.00E+02S 2.00E+02S\n'
        ' H   3 0E0 0E0 1E1 2E1\n SUBTOT 0E0 0E0 1E1 2E1\n'
        ' TOTAL 0E0 0E0 1E1 2E1\n***** JOB FINISHED *****\n')


class PortableTests(unittest.TestCase):
    def test_multiple_zones_remain_visible(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder)/'fort.6'
            source.write_text(BODY)
            result = read_fort6(source)
            import copy
            second = copy.deepcopy(result.tables[0])
            second.zone, second.volume_cm3 = 2, 20
            result.tables.append(second)
            target = Path(folder)/'multi.xlsx'
            export_excel(result, target)
            book = openpyxl.load_workbook(target)
            sheet = book['Activity']
            self.assertEqual(sheet['A7'].value, 'Zone')
            self.assertEqual(sheet.freeze_panes, 'B10')
            self.assertEqual([sheet.cell(7, c).value for c in range(2, 6)], [1, 2, 1, 2])
            book.close()

    def test_real_five_family_values_match_authoritative_payload(self):
        source = Path(__file__).resolve().parents[1]/'6000/fort.6'
        if not source.is_file():
            self.skipTest('Real data are not included in the distribution.')
        result = read_fort6(source)
        payload = workbook_payload(result)
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder)/'6000.xlsx'
            export_excel(result, target)
            book = openpyxl.load_workbook(target, read_only=True)
            for spec in payload['sheets']:
                expected = visible_matrix(payload, spec)
                actual = list(book[spec['name']].iter_rows(min_row=6, values_only=True))
                self.assertEqual(len(actual), len(expected))
                for actual_row, expected_row in zip(actual, expected):
                    self.assertEqual(len(actual_row), len(expected_row))
                    for actual_value, expected_value in zip(actual_row, expected_row):
                        if isinstance(expected_value, (int, float)):
                            import math
                            self.assertTrue(math.isclose(actual_value, expected_value, rel_tol=3e-15, abs_tol=0))
                        else:
                            self.assertEqual(actual_value, expected_value)
            book.close()

    def test_export_view_numeric_values_and_overwrite(self):
        with tempfile.TemporaryDirectory(prefix='project with spaces ') as folder:
            root = Path(folder)
            cell = root/'1234'
            cell.mkdir()
            (cell/'fort.6').write_text(BODY)
            result = read_fort6(cell/'fort.6')
            payload = workbook_payload(result)
            destination = root/'results/1234.xlsx'
            export_excel(result, destination)
            book = openpyxl.load_workbook(destination)
            for spec in payload['sheets']:
                matrix = visible_matrix(payload, spec)
                sheet = book[spec['name']]
                if matrix:
                    actual = list(sheet.iter_rows(min_row=6, values_only=True))
                    self.assertEqual(actual, [tuple(row) for row in matrix])
                    self.assertNotIn('Zone', [row[0] for row in actual])
                    self.assertNotIn('SUBTOT', [row[0] for row in actual])
                    self.assertEqual(sheet.freeze_panes, 'B9')
            self.assertEqual(book['Activity']['B9'].value, 1)
            self.assertEqual(book['Activity']['B9'].number_format, '0.000E+00')
            book.close()
            with self.assertRaises(FileExistsError):
                export_excel(result, destination)
            export_excel(result, destination, overwrite=True)

    def test_change_project_export_zip_and_prevent_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            first, second = root/'first', root/'second'
            for project in (first, second):
                (project/'1234').mkdir(parents=True)
                (project/'1234/fort.6').write_text(BODY)
            store = ResultStore(first)
            store.discover()
            store.load('1234')
            store.change_directory(second)
            self.assertEqual(store.root, second)
            self.assertFalse(store._cache)
            actions = DesktopActions(store)
            response = actions.export(['1234'], False)
            self.assertEqual(response['files'], ['1234.xlsx'])
            self.assertTrue((second/'results/1234.xlsx').is_file())
            self.assertFalse((first/'results').exists())
            import io
            token = response['download'].split('/')[-1]
            with ZipFile(io.BytesIO(actions.downloads[token])) as archive:
                self.assertEqual(archive.namelist(), ['1234.xlsx'])
            with self.assertRaisesRegex(ValueError, 'Already exist'):
                actions.export(['1234'], False)
            with self.assertRaises(ValueError):
                store.change_directory(root/'nonexistent')
            self.assertEqual(store.root, second)

    def test_connections_close_after_grace_but_reload_survives(self):
        sessions = BrowserSessions()
        sessions.update('first', 0, False, now=0)
        sessions.update('second', 0, False, now=0)
        sessions.update('first', 1, True, now=1)
        self.assertFalse(sessions.should_stop(now=10))
        sessions.update('second', 1, True, now=11)
        self.assertFalse(sessions.should_stop(now=11))
        self.assertTrue(sessions.should_stop(now=15))


if __name__ == '__main__':
    unittest.main()

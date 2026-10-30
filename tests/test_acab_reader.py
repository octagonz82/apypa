import math
from pathlib import Path
import tempfile
import unittest

from acab_reader import (AcabResult, ParseError, ResultTable, TimePoint,
                         number, numbers, read_fort6, workbook_payload, time_text)


ROOT = Path(__file__).resolve().parents[1]


class NumberTests(unittest.TestCase):
    def test_fortran_numbers(self):
        self.assertEqual(number('1.25D+03'), 1250)
        self.assertEqual(number('1.25-100'), 1.25e-100)
        self.assertEqual(numbers('1.000E+00-2.000E-02', 2), [1, -0.02])

    def test_invalid_values_are_never_zero(self):
        for value in ('', '*****', 'NaN', '1E999', '1E-999'):
            with self.subTest(value=value), self.assertRaises(ParseError):
                number(value)
        with self.assertRaises(ParseError):
            numbers('0.0 1.0 garbage', 2)
        with self.assertRaises(ParseError):
            numbers('0.0', 2)


class TimeTextTests(unittest.TestCase):
    def test_english_units_and_omission_of_zero_components(self):
        examples = {0: '', 1: '1 second', 107: '107 seconds', 3600: '1 hour',
                    3601: '1 hour and 1 second', 86400: '1 day',
                    2592000: '1 month', 31536000: '1 year', 0.125: '0.125 seconds'}
        for seconds, expected in examples.items():
            with self.subTest(seconds=seconds):
                self.assertEqual(time_text(seconds), expected)
        value = 31536000 + 2*2592000 + 3*86400 + 4*3600 + 5
        self.assertEqual(time_text(value), '1 year, 2 months, 3 days, 4 hours and 5 seconds')

    def test_invalid_durations(self):
        for value in (-1, float('inf'), float('nan')):
            with self.subTest(value=value), self.assertRaises(ValueError):
                time_text(value)


class MatrixTests(unittest.TestCase):
    def result(self, tables, boundaries=None):
        return AcabResult(ROOT/'6000/fort.6', '6000', tables, boundaries or [], [], True)

    def table(self, time_set, rows, times=None):
        return ResultTable('activity', time_set, 1, 10, 1,
                           times or [TimePoint('RESTART', None), TimePoint('1.00E+02S', 100)], rows)

    def test_repeated_times_are_preserved_and_only_all_zero_rows_removed(self):
        result = self.result([
            self.table(6, {'H3': [0, 0], 'Co60': [0, 0], 'Be7': [0, 0], 'TOTAL': [0, 0]}),
            self.table(7, {'H3': [0, 1e-90], 'Co60': [0, 0], 'TOTAL': [0, 1e-90]}),
        ])
        sheet = workbook_payload(result)['sheets'][1]
        self.assertEqual(sheet['retained'], 1)
        self.assertEqual(sheet['omitted'], 2)
        self.assertEqual(sheet['matrix'][0], ['Time set', 6, 7])
        self.assertEqual(sheet['matrix'][2], ['Isotope / time [s]', 100, 100])
        self.assertEqual(sheet['matrix'][3], ['Equivalent duration', '100 seconds', '100 seconds'])
        self.assertEqual(sheet['matrix'][4], ['H3', 0, 1e-90/10])

    def test_absent_observations_stay_blank(self):
        result = self.result([self.table(6, {'H3': [1, 2], 'TOTAL': [1, 2]}),
                              self.table(7, {'Co60': [3, 4], 'TOTAL': [3, 4]})])
        matrix = workbook_payload(result)['sheets'][1]['matrix']
        self.assertEqual(matrix[4], ['H3', 2/10, None])
        self.assertEqual(matrix[5], ['Co60', None, 4/10])

    def test_visible_times_sorted_and_zero_filter_uses_visible_samples(self):
        result = self.result([self.table(7, {'H3': [5, 0, 2], 'TOTAL': [5, 0, 2]},
                                        [TimePoint('INITIAL', None), TimePoint('300S', 300), TimePoint('100S', 100)]),
                              self.table(6, {'Co60': [4, 0], 'TOTAL': [4, 0]})])
        matrix = workbook_payload(result)['sheets'][1]['matrix']
        self.assertEqual(matrix[2], ['Isotope / time [s]', 100, 100, 300])
        self.assertEqual(matrix[4], ['H3', 2/10, None, 0])
        self.assertNotIn('Co60', [row[0] for row in matrix[4:]])

    def test_gamma_orientation_and_printed_totals(self):
        table = ResultTable('gamma', 6, 1, 10, 1,
            [TimePoint('RESTART', None), TimePoint('1E2S', 100)],
            {'1E0': [1, 2], '2.5E-1': [3, 4], 'TOTAL': [4.01, 6.01], 'MEV/SEC': [2.5, 4]}, [1, 0.25])
        sheet = workbook_payload(self.result([table], [1.5, 0.5, 0]))['sheets'][4]
        matrix = sheet['matrix']
        self.assertEqual(matrix[0][4:6], [0.5, 1.5])
        self.assertEqual(len(matrix), 2)
        self.assertEqual(matrix[1], [6, 1, 100, '100 seconds', 4, 2, 6.01, 4])
        self.assertEqual(sheet['unit'], 'photons/(cm³·s)')
        self.assertEqual(table.rows['1E0'], [1, 2])  # Raw data are never modified.

    def test_normalization_uses_each_zone_volume_and_totals(self):
        first = self.table(6, {'H3': [0, 100], 'SUBTOT': [0, 100], 'TOTAL': [0, 120]})
        second = self.table(7, {'H3': [0, 100], 'SUBTOT': [0, 100], 'TOTAL': [0, 120]})
        second.zone, second.volume_cm3 = 2, 20
        sheet = workbook_payload(self.result([first, second]))['sheets'][1]
        self.assertEqual(sheet['unit'], 'Bq/cm³')
        self.assertEqual(sheet['matrix'][4:], [['H3', 10, 5], ['TOTAL', 12, 6]])
        self.assertIn('SUBTOT', first.rows)

    def test_physical_units_and_no_volume_division_for_dose(self):
        for family, unit, expected in [('composition', 'mol/cm³', 10),
                                       ('heat', 'W/cm³', 10), ('dose', 'Sv/h', 100)]:
            table = self.table(6, {'H3': [0, 100], 'TOTAL': [0, 120]})
            table.family = family
            sheet = next(s for s in workbook_payload(self.result([table]))['sheets']
                         if s['family'] == family)
            self.assertEqual(sheet['unit'], unit)
            self.assertEqual(sheet['matrix'][4], ['H3', expected])
            self.assertEqual(table.rows['H3'], [0, 100])

    def test_invalid_volume_and_missing_gamma_limits_fail(self):
        table = self.table(6, {'H3': [0, 1], 'TOTAL': [0, 1]})
        for volume in (None, 0, -1, float('nan')):
            table.volume_cm3 = volume
            with self.subTest(volume=volume), self.assertRaises(ParseError):
                workbook_payload(self.result([table]))
        table.family, table.volume_cm3 = 'gamma', 10
        with self.assertRaisesRegex(ParseError, 'gamma boundaries'):
            workbook_payload(self.result([table]))


class ReaderTests(unittest.TestCase):
    def fixture(self, body):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'fort.6'
            path.write_text(body, encoding='utf8')
            return read_fort6(path)

    HEADER = ('1Test 6. TIME SET\n VOLUME OF ZONE : 10 CCM\n'
              '0 NUCLIDE RADIOACTIVITY, DISINTEGRATIONS/SEC\n ZONE 1\n\n'
              ' INITIAL RESTART 1.00E+02S\n')

    def test_states_and_meta_stable_labels(self):
        result = self.fixture(self.HEADER + ' CO 60M1.000D+00 2.000E+00 3.000E+00\n'
                              ' TOTAL 1.000E+00 2.000E+00 3.000E+00\n***** JOB FINISHED *****\n')
        self.assertEqual(result.tables[0].rows['Co60m'], [1, 2, 3])
        self.assertEqual([t.seconds for t in result.tables[0].times], [None, None, 100])
        self.assertTrue(result.finished)

    def test_truncated_table_fails(self):
        with self.assertRaises(ParseError):
            self.fixture(self.HEADER + ' H   3 1E0 2E0 3E0\n')

    def test_invalid_cell_reports_source_line(self):
        with self.assertRaisesRegex(ParseError, r'fort\.6:7:'):
            self.fixture(self.HEADER + ' H   3 1E0 ***** 3E0\n TOTAL 1E0 2E0 3E0\n')

    def test_all_ten_real_files_against_printed_rows(self):
        paths = sorted(ROOT.glob('*/fort.6'))
        if not paths:
            self.skipTest("The user's real data are not included in the repository.")
        self.assertEqual(len(paths), 10)
        for path in paths:
            with self.subTest(cell=path.parent.name):
                result = read_fort6(path)
                source = path.read_text().splitlines()
                self.assertTrue(result.finished)
                self.assertEqual(len(result.tables), 30)
                self.assertEqual(len(result.gamma_boundaries_mev), 25)
                for family in ('composition', 'activity', 'heat', 'dose', 'gamma'):
                    tables = result.family(family)
                    self.assertEqual([t.time_set for t in tables], [6, 7, 8, 9, 10, 11])
                    for table in tables:
                        # Independent float/split decoding of every printed data row.
                        seen = {}
                        for line in source[table.line:]:
                            if family == 'gamma':
                                parts = line.split()
                                if parts and parts[0] == '0':
                                    parts = parts[1:]
                                if not parts or parts[0] in ('(MEV)', 'EMEAN', 'ZONE'):
                                    continue
                                if parts[0] not in table.rows:
                                    continue
                                key, values = parts[0], [float(x) for x in parts[1:]]
                            else:
                                key = line[:7].replace(' ', '')
                                if key not in ('TOTAL', 'SUBTOT'):
                                    key = key.capitalize()
                                if key not in table.rows:
                                    continue
                                values = [float(x) for x in line[7:].split()]
                            self.assertEqual(values, table.rows[key])
                            seen[key] = values
                            if key == ('MEV/SEC' if family == 'gamma' else 'TOTAL'):
                                break
                        self.assertEqual(set(seen), set(table.rows))
                        # Printed totals may differ from sums due to ACAB rounding.
                        detail = [row for k, row in table.rows.items() if k not in ('TOTAL', 'SUBTOT', 'MEV/SEC')]
                        reference = table.rows.get('SUBTOT', table.rows['TOTAL'])
                        tolerance = 0.015 if family == 'gamma' else 0.002
                        for j, total in enumerate(reference):
                            self.assertTrue(math.isclose(math.fsum(row[j] for row in detail), total,
                                                         rel_tol=tolerance, abs_tol=1e-95))


if __name__ == '__main__':
    unittest.main()

import math
from pathlib import Path
import tempfile
import unittest
import io
import json
from contextlib import redirect_stdout
from unittest.mock import patch

from acab_reader import read_fort6
from apypa_plotter import plot_data, ResultStore, existing_viewer, main, BrowserSessions

ROOT = Path(__file__).resolve().parents[1]


class PlotterTests(unittest.TestCase):
    def test_last_page_close_and_reload_grace(self):
        sessions = BrowserSessions(grace=3, timeout=300)
        self.assertFalse(sessions.should_stop(now=1000))  # Wait for first page.
        sessions.update('a', 1, False, now=1000)
        sessions.update('b', 1, False, now=1000)
        sessions.update('a', 2, True, now=1001)
        self.assertFalse(sessions.should_stop(now=1005))  # Other tab still open.
        sessions.update('b', 2, True, now=1006)
        self.assertFalse(sessions.should_stop(now=1006))
        sessions.update('reload', 1, False, now=1007)
        self.assertFalse(sessions.should_stop(now=1010))
        sessions.update('reload', 3, True, now=1011)
        sessions.update('reload', 2, False, now=1012)  # Late heartbeat must not reopen.
        self.assertFalse(sessions.should_stop(now=1012))
        self.assertTrue(sessions.should_stop(now=1015))

    def test_abrupt_browser_exit_expires(self):
        sessions = BrowserSessions(grace=3, timeout=300)
        sessions.update('page', 1, False, now=0)
        self.assertFalse(sessions.should_stop(now=299))
        self.assertFalse(sessions.should_stop(now=300))
        self.assertTrue(sessions.should_stop(now=303))

    def test_reuse_only_same_app_and_directory(self):
        identity = {'app': 'apypa-plotter', 'root': str(ROOT.resolve())}
        for payload, expected in [(identity, identity),
                                  ({**identity, 'root': '/other'}, None),
                                  ({'cells': []}, None)]:
            with patch('apypa_plotter.urllib.request.urlopen',
                       return_value=io.BytesIO(json.dumps(payload).encode())):
                self.assertEqual(existing_viewer(8765, ROOT), expected)

    def test_repeated_start_does_not_spawn_server(self):
        with patch('apypa_plotter.existing_viewer', return_value={'app': 'apypa-plotter'}), \
             patch('apypa_plotter.subprocess.Popen') as spawn, \
             patch('apypa_plotter.ThreadingHTTPServer') as server, \
             redirect_stdout(io.StringIO()) as output:
            self.assertEqual(main([str(ROOT), '--background']), 0)
            self.assertIn('http://127.0.0.1:8765', output.getvalue())
            spawn.assert_not_called()
            server.assert_not_called()

    @classmethod
    def setUpClass(cls):
        if not (ROOT/'6000/fort.6').is_file():
            raise unittest.SkipTest("The user's real data are not included in the repository.")
        cls.result = read_fort6(ROOT/'6000/fort.6')
        cls.data = plot_data(cls.result)

    def test_physical_units_and_real_data(self):
        f = self.data['families']
        self.assertEqual([f[k]['unit'] for k in ('composition','activity','heat','dose','gamma')],
                         ['mol/cm³','Bq/cm³','W/cm³','Sv/h','photons/(cm³·s)'])
        activity = next(s['values'] for s in f['activity']['series'] if s['name']=='H3')
        self.assertTrue(math.isclose(activity[0],38470000/610712,rel_tol=1e-12))
        dose = next(s['values'] for s in f['dose']['series'] if s['name']=='Be7')
        self.assertEqual(dose[0],6.363e-11)
        for data in f.values():
            self.assertNotIn('SUBTOT', [s['name'] for s in data['series']])
            self.assertEqual(len(data['samples']),60)
            self.assertTrue(all(isinstance(s['time'],(int,float)) for s in data['samples']))
            self.assertTrue(all(s['label'] for s in data['samples']))

    def test_gamma_boundaries_order_and_no_additional_normalization(self):
        gamma = self.data['families']['gamma']
        self.assertEqual(gamma['bins'][0],{'lower':0,'upper':.01})
        self.assertEqual(gamma['bins'][-1],{'lower':14,'upper':20})
        self.assertEqual(len(gamma['bins']),24)
        self.assertEqual(gamma['values'][0][0],57.4)
        self.assertEqual(gamma['totals'][0],409000)

    def test_cache_invalidation_and_unknown_cells(self):
        # Small isolated copy makes this test independent of user source changes.
        with tempfile.TemporaryDirectory() as directory:
            cell=Path(directory)/'1234';cell.mkdir()
            source=cell/'fort.6'
            body=('1Test 6. TIME SET\n VOLUME OF ZONE : 10 CCM\n'
                  '0 NUCLIDE RADIOACTIVITY, DISINTEGRATIONS/SEC\n ZONE 1\n'
                  ' INITIAL RESTART 1.00E+02S\n H   3 0E0 0E0 1E1\n'
                  ' TOTAL 0E0 0E0 1E1\n***** JOB FINISHED *****\n')
            source.write_text(body)
            store=ResultStore(Path(directory))
            self.assertEqual(store.discover()[0]['id'],'1234')
            first=store.load('1234');self.assertIs(store.load('1234'),first)
            self.assertEqual(first['families']['gamma']['bins'],[])
            source.write_text(body.replace('1E1','2E1')+'\n')
            self.assertIsNot(store.load('1234'),first)
            with self.assertRaises(KeyError):store.load('../1234')


if __name__=='__main__':unittest.main()

"""
Failures inside liborbfit raise BKOrbitError instead of ending the process,
concurrent calls from threads give the same answers as serial ones, and the
library exports only its public entry points.
"""
import os
import tempfile
import threading
import unittest

from astropy.time import Time

import mp_ephem
from mp_ephem import bk_orbit

__PATH__ = os.path.dirname(__file__)


class LibraryExports(unittest.TestCase):

    def test_only_api_exported(self):
        lib = bk_orbit._load_orbfit()
        for name in ('fitradec', 'predict', 'predict_helio', 'abg_to_aei', 'orbfit_last_error',
                     'set_ephem_file', 'set_observatory_file', 'set_mpc_dtheta', 'utc_jd_to_tt'):
            self.assertTrue(hasattr(lib, name), name)
        for name in ('gaussj', 'dvector', 'nrerror', 'kbo3d', 'read_abg', 'orbfit_fail'):
            self.assertFalse(hasattr(lib, name), name)


class LibraryErrors(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        observations = mp_ephem.EphemerisReader().read(os.path.join(__PATH__, 'data', 'o3o08.mpc'))
        cls.orbit = mp_ephem.BKOrbit(observations)

    def test_date_outside_ephemeris(self):
        with self.assertRaisesRegex(mp_ephem.BKOrbitError, 'out of range of ephemeris'):
            self.orbit.predict(Time('2300-01-01', scale='utc'))
        with self.assertRaisesRegex(mp_ephem.BKOrbitError, 'out of range of ephemeris'):
            self.orbit.predict_helio(Time('2300-01-01', scale='utc'))

    def test_recovers_after_failure(self):
        date = Time('2016-04-09', scale='utc')
        self.orbit.predict(date)
        expected = self.orbit.coordinate.ra.degree
        with self.assertRaises(mp_ephem.BKOrbitError):
            self.orbit.predict(Time('2300-01-01', scale='utc'))
        self.orbit.predict(date)
        self.assertEqual(self.orbit.coordinate.ra.degree, expected)

    def test_bad_abg_file(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.abg') as abg:
            abg.write('# not an orbit\n1 2 3\n')
            abg.flush()
            with self.assertRaisesRegex(mp_ephem.BKOrbitError, 'a/b/g'):
                mp_ephem.BKOrbit(None, os.path.join(__PATH__, 'data', 'o3o08.mpc'), abg_file=abg.name)

    def test_threads_match_serial(self):
        dates = [Time(2457000.5 + 40 * i, format='jd', scale='utc') for i in range(8)]
        serial = []
        for date in dates:
            self.orbit.predict(date)
            serial.append((self.orbit.ra.degree, self.orbit.dec.degree))

        results = {}
        errors = []

        def work(index):
            try:
                orbit = mp_ephem.BKOrbit(None, os.path.join(__PATH__, 'data', 'o3o08.mpc'))
                for _ in range(3):
                    orbit.predict(dates[index])
                results[index] = (orbit.ra.degree, orbit.dec.degree)
            except Exception as ex:  # pragma: no cover - reported below
                errors.append(ex)

        threads = [threading.Thread(target=work, args=(i,)) for i in range(len(dates))]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(errors, [])
        for index, (ra, dec) in enumerate(serial):
            self.assertAlmostEqual(results[index][0], ra, delta=1e-9)
            self.assertAlmostEqual(results[index][1], dec, delta=1e-9)

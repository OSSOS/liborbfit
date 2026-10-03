import ctypes
import glob
import os
import unittest

import mp_ephem
from astropy.time import Time


class TimeTest(unittest.TestCase):

    def test_precision(self):
        time_str = '2001 01 01.000010'
        time_precision = 6
        this_time = Time(time_str, scale='utc', precision=time_precision)
        # print(this_time.isot)
        self.assertEqual(str(this_time), time_str)


class UtcToTtTest(unittest.TestCase):
    """DE405 is evaluated at TT. The UTC offset must follow leap seconds."""

    @classmethod
    def setUpClass(cls):
        lib = glob.glob(os.path.join(os.path.dirname(mp_ephem.__file__), 'orbfit*.so'))[0]
        cls.orbfit = ctypes.CDLL(lib)
        cls.orbfit.utc_jd_to_tt.restype = ctypes.c_double
        cls.orbfit.utc_jd_to_tt.argtypes = [ctypes.c_double]

    def assert_matches_astropy(self, utc_string):
        utc = Time(utc_string, scale='utc')
        got = self.orbfit.utc_jd_to_tt(utc.jd)
        # A tenth of a millisecond is far inside the ephemeris time tolerance.
        self.assertLess(abs(got - utc.tt.jd) * 86400.0, 1.0e-4,
                        msg='{}: C TT {} vs astropy {}'.format(utc_string, got, utc.tt.jd))

    def test_leap_second_boundaries(self):
        for stamp in (
            '1960-01-01',
            '1965-06-15',
            '1971-12-31',
            '1972-01-01',
            '1961-07-31 18:00:00',
            '1998-12-31 12:00:00',
            '1998-12-31 23:59:59',
            '1999-01-01',
            '2005-12-31 23:59:59',
            '2006-01-01',
            '2016-12-31 12:00:00',
            '2016-12-31 23:59:59',
            '2016-12-31 23:59:60',
            '2017-01-01',
            '2013-04-03 15:06:08',
            '2026-10-03',
        ):
            self.assert_matches_astropy(stamp)

    def test_modern_offset_is_not_the_1999_constant(self):
        utc = Time('2026-10-03', scale='utc')
        offset_seconds = (self.orbfit.utc_jd_to_tt(utc.jd) - utc.jd) * 86400.0
        self.assertAlmostEqual(offset_seconds, 69.184, places=3)
        self.assertGreater(abs(offset_seconds - 64.0), 5.0)

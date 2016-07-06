import unittest
from mp_ephem import time_mpc
from astropy.time import Time


class TimeTest(unittest.TestCase):

    def test_precision(self):
        time_str = '2001 01 01.00001'
        time_precision = 5
        this_time = Time(time_str, scale='utc', precision=time_precision)
        self.assertEquals(str(this_time), time_str)
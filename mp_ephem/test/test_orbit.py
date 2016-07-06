import unittest
import mp_ephem
from astropy import units
import os


class OrbitFit(unittest.TestCase):

    def setUp(self):
        mpc_filename = 'data/o3o08.mpc'
        self.observations = mp_ephem.EphemerisReader().read(mpc_filename)

    def test_orbit(self):
        """
        Test that the Keplarian orbit elements returned by BJOrbit matches the expected values

        :return:
        """
        self.orbit = mp_ephem.BKOrbit(self.observations)
        self.assertAlmostEqual(self.orbit.a, 33.3 * units.AU)
        self.assertAlmostEqual(self.orbit.e, 0.1)
        self.assertAlmostEqual(self.orbit.inc, 1.00 * units.degree)
        self.assertAlmostEqual(self.orbit.Node, 180.0 * units.degree)
        self.assertAlmostEqual(self.orbit.om, 180.0 * units.degree)
        self.assertAlmostEqual(self.orbit.T, 2450000.0 * units.day)
        self.assertAlmostEqual(self.orbit.epoch, 2450000.0 * units.day)

    def test_data(self):
        """
        Test that Ephemeris Binary and observtories file form JPL can be found by the environment variable.
        :return:
        """

        print os.environ['ORBIT_EPHEMERIS']
        self.assertTrue(os.access(os.environ['ORBIT_EPHEMERIS'], os.R_OK))
        self.assertTrue(os.access(os.environ['ORBIT_OBSERVATORIES'], os.R_OK))

    def test_orbfit_residuals(self):
        mpc_lines = ("     HL7j2    C2013 04 03.62926 17 12 01.16 +04 13 33.3          24.1 R      568",
                     "     HL7j2    C2013 04 04.58296 17 11 59.80 +04 14 05.5          24.0 R      568",
                     "     HL7j2    C2013 05 03.52252 17 10 38.28 +04 28 00.9          23.4 R      568",
                     "     HL7j2    C2013 05 08.56725 17 10 17.39 +04 29 47.8          23.4 R      568")

        observations = []
        for line in mpc_lines:
            observations.append(mp_ephem.Observation.from_string(line))

        this_orbit = mp_ephem.BKOrbit(observations=observations)

        for observation in observations:
            this_orbit.predict(observation.date, 568)
            self.assertLess(observation.ra_residual, 0.3)
            self.assertLess(observation.dec_residual, 0.3)

    def test_null_obseravtion(self):
        mpc_lines = ("!    HL7j2    C2013 04 03.62926 17 12 01.16 +04 13 33.3          24.1 R      568",
                     "     HL7j2    C2013 04 04.58296 17 11 59.80 +04 14 05.5          24.0 R      568",
                     "     HL7j2    C2013 05 03.52252 17 10 38.28 +04 28 00.9          23.4 R      568",
                     "     HL7j2    C2013 05 08.56725 17 10 17.39 +04 29 47.8          23.4 R      568")

        observations = []
        for line in mpc_lines:
            observations.append(mp_ephem.Observation.from_string(line))

        this_orbit = mp_ephem.BKOrbit(observations)
        self.assertAlmostEqual(this_orbit.a, 135.75, 1)

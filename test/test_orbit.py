import unittest
import mp_ephem
from astropy import units
import os


class OrbitFit(unittest.TestCase):

    def setUp(self):
        mpc_filename = 'data/o3o08.mpc'
        self.abg_filename = 'data/o3o08.abg'
        self.observations = mp_ephem.EphemerisReader().read(mpc_filename)
        self.orbit = mp_ephem.BKOrbit(self.observations)

    def test_orbit(self):
        """
        Test that the Keplarian orbit elements returned by BJOrbit matches the expected values

        :return:
        """
        self.assertAlmostEqual(self.orbit.a.to(units.AU).value, 39.3419, 3)
        self.assertAlmostEqual(self.orbit.e.value, 0.2778, 3)
        self.assertAlmostEqual(self.orbit.inc.to(units.degree).value, 8.05, 2)
        self.assertAlmostEqual(self.orbit.Node.to(units.degree).value, 113.85, 2)
        self.assertAlmostEqual(self.orbit.om.to(units.degree).value, 66.24, 2)
        self.assertAlmostEqual(self.orbit.T.to(units.day).value, 2447884.5070, 3)
        self.assertAlmostEqual(self.orbit.epoch.jd, 2456392.05115, 4)

    def test_data(self):
        """
        Test that Ephemeris Binary and observtories file form JPL can be found by the environment variable.
        :return:
        """

        self.assertTrue(os.access(os.environ['ORBIT_EPHEMERIS'], os.R_OK))
        self.assertTrue(os.access(os.environ['ORBIT_OBSERVATORIES'], os.R_OK))

    def test_abg_load(self):
        """
        Test that loading an abg file returns the same results as calling fit_radec
        :return:
        """

        orbit1 = mp_ephem.BKOrbit(self.observations, abg_file=self.abg_filename)
        orbit2 = mp_ephem.BKOrbit(self.observations)
        for attr in ['a', 'e', 'Node', 'inc', 'om', 'T', 'distance']:
            self.assertEqual(getattr(orbit1, attr),
                             getattr(orbit2, attr))


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
            this_orbit.compute_residuals()
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
        self.assertAlmostEqual(this_orbit.a.to(units.au).value, 137.91, 1)

# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import unittest
import mp_ephem
from astropy import units, coordinates
import os
from mp_ephem.ephem import ObserverLocation

__PATH__ = os.path.dirname(__file__)


class HSTFormat(unittest.TestCase):

    def setUp(self):
        self.mpc_filename = os.path.join(__PATH__, 'data/hst.mpc')

    def test_parser(self):
        obs = mp_ephem.EphemerisReader().read(self.mpc_filename)
        self.assertIsInstance(obs[0].location, ObserverLocation)

    def test_fitradec(self):
        orbit = mp_ephem.BKOrbit(None, self.mpc_filename)
        orbit.predict(orbit.observations[0].date)
        self.assertAlmostEqual(orbit.a.to('au').value, 17.08, 0)

class SimonFormat(unittest.TestCase):

    def setUp(self):
        self.mpc_filename = os.path.join(__PATH__,'data/simon_format.txt')

    def test_parser(self):
        obs = mp_ephem.EphemerisReader().read(self.mpc_filename)
        self.assertIsInstance(obs[0], mp_ephem.Observation)

    def test_fitradec(self):
        orbit = mp_ephem.BKOrbit(None, self.mpc_filename)
        orbit.predict(orbit.observations[0].date)
        self.assertAlmostEqual(orbit.a.to('au').value, 43.8015219, 3)


class OrbitFit(unittest.TestCase):

    def setUp(self):
        mpc_filename = os.path.join(__PATH__,'data/o3o08.mpc')
        self.abg_filename = os.path.join(__PATH__, 'data/o3o08.abg')
        self.observations = mp_ephem.EphemerisReader().read(mpc_filename)
        self.orbit = mp_ephem.BKOrbit(self.observations)
        self.mpc_lines = ("     HL7j2    C2013 04 03.62926 17 12 01.16 +04 13 33.3          24.1 R      568",
                          "     HL7j2    C2013 04 04.58296 17 11 59.80 +04 14 05.5          24.0 R      568",
                          "     HL7j2    C2013 05 03.52252 17 10 38.28 +04 28 00.9          23.4 R      568",
                          "     HL7j2    C2013 05 08.56725 17 10 17.39 +04 29 47.8          23.4 R      568")
        observations = []
        for line in self.mpc_lines:
            observations.append(mp_ephem.ObsRecord.from_string(line))

        self.example_orbit = mp_ephem.BKOrbit(observations=observations)


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
        self.assertAlmostEqual(self.orbit.T.to(units.day).value, 2447884.5762, 3)
        self.assertAlmostEqual(self.orbit.epoch.jd, 2456392.05115, 4)

    def test_summarize(self):
        """
        Print a summary of the observation.
        :return:
        """
        self.assertIsInstance(self.orbit.summarize(), str)

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

    def test_OSSOSParser(self):
        mpc_line = " O13BL3UV     C2013 08 02.50855 01 00 04.549+04 59 01.53         24.3 r      568 O 1645236p27 L3UV Y 106.35 4301.85 0.20 0 24.31 0.15 % hurrah!"
        obs = mp_ephem.Observation.from_string(mpc_line)
        self.assertIsInstance(obs.comment, mp_ephem.OSSOSComment)
        self.assertAlmostEqual(obs.comment.x, 106.35)
        self.assertAlmostEqual(obs.comment.y, 4301.85)
        mpc_line = "     K01QX1F 1C2000 08 26.21908 23 10 50.37 -05 45 16.1          22.6 R      807 19000101_500_1 20170607 0000000000                       link o3l06PD = 2001 QF331 = K01QX1F"
        obs = mp_ephem.Observation.from_string(mpc_line)
        self.assertIsInstance(obs.comment, mp_ephem.MPCComment)
        mpc_line = " O13BL3SX     C2013 09 29.42805 00 46 59.117+02 09 12.81         24.1 r      568 20130929_568_1 20180216 1000000000                      O   1656905p14 L3SX        Y  1873.54 4214.87 0.11 2 24.07 0.14 %"
        mpc_line = " O13BL3SX    EC2013 12 05.23284 00 42 51.561+01 44 04.10                     568 20130802_568_1 20140817 0000000000                      O 1672595p13 O13BL3SX ZE  170.4 4611.6            UUUU % just visible"
        obs = mp_ephem.Observation.from_string(mpc_line)
        self.assertIsInstance(obs.comment, mp_ephem.OSSOSComment)
        mpc_line = " O13BL3SX     C2014 06 25.58497 00 55 37.028+03 03 50.66         23.6 r      568 20140102_568_1 20140817 0000000000                      O 1722362p11 O13BL3SX Y  1905.8 1009.5 23.57 0.10 UUUU % "
        obs = mp_ephem.Observation.from_string(mpc_line)
        self.assertIsInstance(obs.comment, mp_ephem.OSSOSComment)
        mpc_line = "q3615K06UW1O HC2018 09 12.62432 01 13 26.144+05 53 05.74               q~2kJo568"
        obs = mp_ephem.Observation.from_string(mpc_line)
        self.assertIsInstance(obs.comment, str)

    def dont_test_orbfit_residuals(self):
        for observation in self.observations:
            self.example_orbit.predict(observation.date, 568)
            self.example_orbit.compute_residuals()
            self.assertLess(self.example_orbit.observations[0].ra_residual, 0.3)
            self.assertLess(self.example_orbit.observations[0].dec_residual, 0.3)

    def test_predict_helio(self):
        """
        Ensure that the geocentric and heliocentric coordinates transform correctly.
        :return:
        """
        self.orbit.predict("2013 04 09.55115")
        #print(self.orbit.heliocentric)
        #print(self.orbit.geocentric)
        #print(self.orbit.geocentric.transform_to('geocentricmeanecliptic'))
        #print(self.orbit.coordinate)

    def test_null_obseravtion(self):
        self.assertAlmostEqual(self.example_orbit.a.to(units.au).value/100, 137.91/100, 1)

    def test_alphanumeric_obscode(self):
        """
        T14 (CFHT) sits beside 568 on Maunakea, so both codes must give the same orbit.
        """
        observations = [obs for obs in self.observations if not obs.null_observation]
        for obs in observations:
            obs.observatory_code = 'T14'
        self.assertEqual(observations[0].observatory_code, 'T14')
        orbit = mp_ephem.BKOrbit(observations)
        self.assertAlmostEqual(orbit.a.to(units.AU).value, self.orbit.a.to(units.AU).value, 3)
        self.assertAlmostEqual(orbit.e.value, self.orbit.e.value, 4)

    def test_predict_space_observatory(self):
        """
        A space observatory without a supplied position is predicted from the geocenter.
        """
        self.orbit.predict(self.observations[0].date, obs_code=250)
        geocentric = self.orbit.coordinate
        self.orbit.predict(self.observations[0].date, obs_code=500)
        self.assertAlmostEqual(geocentric.separation(self.orbit.coordinate).to(units.arcsec).value, 0, 6)

    def test_tnodb_discovery_flags(self):
        orbit = mp_ephem.BKOrbit(None, ast_filename=os.path.join(__PATH__,'data/o4h29.ast'))
        for observation in orbit.observations:
            self.assertTrue(observation.discovery)

class ObservatoryCodes(unittest.TestCase):

    def test_obscode_to_int(self):
        for code, value in ((568, 568), ('568', 568), ('000', 0), ('I11', 1811), ('T14', 2914),
                            ('W84', 3284), ('Z99', 3599), ('a05', 3605)):
            self.assertEqual(mp_ephem.ephem.obscode_to_int(code), value)
        for code in ('T1', 'TT4', ''):
            self.assertRaises(ValueError, mp_ephem.ephem.obscode_to_int, code)

    def test_known_codes(self):
        for code in ('500', '568', 'T14', 'W84', '250'):
            self.assertEqual(mp_ephem.ObsRecord(observatory_code=code).observatory_code, code)


class CLASSYFORM(unittest.TestCase):

    def setUp(self):
        self.mpc_filename = os.path.join(__PATH__, 'data/classy.mpc')

    def test_parser(self):
        obs = mp_ephem.EphemerisReader().read(self.mpc_filename)
        self.assertIsInstance(obs[0], mp_ephem.Observation)
        self.assertIsInstance(obs[0].comment, mp_ephem.OSSOSComment)
        self.assertAlmostEqual(obs[0].comment.likelihood, 12)

    def test_fitradec(self):
        orbit = mp_ephem.BKOrbit(None, self.mpc_filename)
        orbit.predict(orbit.observations[0].date)
        self.assertAlmostEqual(orbit.a.to('au').value, 44.48, 2)  # Example value, adjust as needed
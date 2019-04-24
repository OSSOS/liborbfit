# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import unittest
import mp_ephem
from astropy import units
import os


class SimonFormat(unittest.TestCase):

    def setUp(self):
        self.mpc_filename = 'data/simon_format.txt'

    def test_parser(self):
        obs = mp_ephem.EphemerisReader().read(self.mpc_filename)
        self.assertIsInstance(obs[0], mp_ephem.Observation)

    def test_fitradec(self):
        orbit = mp_ephem.BKOrbit(None, self.mpc_filename)
        orbit.predict(orbit.observations[0].date)
        self.assertAlmostEqual(orbit.a.to('au').value, 43.8026798, 5)

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

    def test_orbfit_residuals(self):
        mpc_lines = ("     HL7j2    C2013 04 03.62926 17 12 01.16 +04 13 33.3          24.1 R      568",
                     "     HL7j2    C2013 04 04.58296 17 11 59.80 +04 14 05.5          24.0 R      568",
                     "     HL7j2    C2013 05 03.52252 17 10 38.28 +04 28 00.9          23.4 R      568",
                     "     HL7j2    C2013 05 08.56725 17 10 17.39 +04 29 47.8          23.4 R      568")

        observations = []
        for line in mpc_lines:
            observations.append(mp_ephem.ObsRecord.from_string(line))

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
            observations.append(mp_ephem.ObsRecord.from_string(line))

        this_orbit = mp_ephem.BKOrbit(observations)
        self.assertAlmostEqual(this_orbit.a.to(units.au).value, 137.91, 1)

    def test_tnodb_discovery_flags(self):
        orbit = mp_ephem.BKOrbit(None, ast_filename='data/o4h29.ast')
        for observation in orbit.observations:
            self.assertTrue(observation.discovery)

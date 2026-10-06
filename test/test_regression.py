"""
Pin the orbit fit, prediction and uncertainty for o3o08 so that changes to the
C library that move results are noticed.

Tolerances allow for input precision at the 0.05 arcsec level (the size of the
float-to-double fix in RA/Dec parsing) but catch any real change in the fit.
"""
import os
import unittest

from astropy import units

import mp_ephem

__PATH__ = os.path.dirname(__file__)


class O3o08Regression(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        observations = mp_ephem.EphemerisReader().read(os.path.join(__PATH__, 'data', 'o3o08.mpc'))
        cls.orbit = mp_ephem.BKOrbit(observations)

    def test_elements(self):
        orbit = self.orbit
        self.assertAlmostEqual(orbit.a.to(units.au).value, 39.34188, delta=2e-3)
        self.assertAlmostEqual(orbit.e.value, 0.277796, delta=1e-4)
        self.assertAlmostEqual(orbit.inc.to(units.degree).value, 8.048165, delta=1e-3)
        self.assertAlmostEqual(orbit.Node.to(units.degree).value, 113.85162, delta=2e-3)
        self.assertAlmostEqual(orbit.om.to(units.degree).value, 66.2447, delta=0.05)
        self.assertAlmostEqual(orbit.T.to(units.day).value, 2447884.58, delta=5)
        self.assertAlmostEqual(orbit.epoch.jd, 2456392.05115, delta=1e-6)

    def test_distance_and_uncertainty(self):
        orbit = self.orbit
        self.assertAlmostEqual(orbit.distance.to(units.au).value, 31.65719, delta=1e-3)
        self.assertAlmostEqual(orbit.da.to(units.au).value, 0.02748, delta=0.003)

    def test_predict_at_first_observation(self):
        orbit = self.orbit
        orbit.predict(orbit.observations[0].date)
        self.assertAlmostEqual(orbit.coordinate.ra.degree, 238.9183567, delta=0.1 / 3600)
        self.assertAlmostEqual(orbit.coordinate.dec.degree, -13.4157912, delta=0.1 / 3600)
        self.assertAlmostEqual(orbit.dra.to(units.arcsec).value, 0.0994, delta=0.01)
        self.assertAlmostEqual(orbit.ddec.to(units.arcsec).value, 0.0630, delta=0.01)

    def test_predict_two_years_out(self):
        orbit = self.orbit
        orbit.predict("2016 04 09.55115")
        self.assertAlmostEqual(orbit.coordinate.ra.degree, 245.3271649, delta=2.0 / 3600)
        self.assertAlmostEqual(orbit.coordinate.dec.degree, -15.1861396, delta=2.0 / 3600)
        self.assertAlmostEqual(orbit.dra.to(units.arcsec).value, 2.382, delta=0.2)


if __name__ == '__main__':
    unittest.main()

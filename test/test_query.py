from unittest import TestCase
from mp_ephem import horizons


class TestQuery(TestCase):

    def test_target(self):
        body = horizons.Body('2002 MS4')
        self.assertIsInstance(body, horizons.Body)
        mag = body.mag
        self.fail()

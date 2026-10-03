from unittest import TestCase

from astropy import units
from astropy.time import Time

from mp_ephem import horizons

# Minimal Horizons CSV response captured for 2002 MS4 on 2020-01-01/02.
# Keeps this test independent of live JPL queries and Time.now().
HORIZONS_FIXTURE = [
    b"****",
    b"   A= 41.96777612227916",
    b" Date__(UT)__HR:MN:SS, Date_________JDUT, , , R.A._(ICRF), DEC_(ICRF),"
    b"  dRA*cosD, d(DEC)/dt,  RA_3sigma, DEC_3sigma,    APmag,  S-brt,"
    b"  SMAA_3sig, SMIA_3sig,   Theta, Area_3sig,                r,       rdot,"
    b"             delta,     deldot,       phi,  PAB-LON,  PAB-LAT,     S-O-T,/r,",
    b"***",
    b"$$SOE",
    b" 2020-Jan-01 00:00:00, 2458849.500000000,*,m,   281.94116,   -6.97964,"
    b"  3.531514,  0.329693,      0.005,      0.003,   20.441,   n.a.,"
    b"    0.00502,   0.00261,  -7.980, 0.0000413,   46.55184511243, -0.3665496,"
    b"  47.4954884093131,  1.0730811,    0.3365, 282.3613,  16.1159,   16.1532,/T,",
    b" 2020-Jan-02 00:00:00, 2458850.500000000,*,m,   281.96516,   -6.97745,"
    b"  3.532164,  0.345306,      0.005,      0.003,   20.440,   n.a.,"
    b"    0.00506,   0.00261,  -7.947, 0.0000416,   46.55163339978, -0.3665946,"
    b"  47.4958693813771,  0.5839054,    0.3339, 282.3753,  16.1162,   16.0271,/T,",
    b"$$EOE",
    b"# obs: 496 (1954-2023)",
]


class TestQuery(TestCase):

    def test_target(self):
        start = Time("2020-01-01", scale="utc")
        stop = start + 1 * units.day
        body = horizons.Body("2002 MS4", start_time=start, stop_time=stop)
        body._data = HORIZONS_FIXTURE
        body.predict(start)

        self.assertIsInstance(body, horizons.Body)
        self.assertAlmostEqual(body.mag, 20.441, places=2)

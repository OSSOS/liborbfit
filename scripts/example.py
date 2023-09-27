from mp_ephem import *
from astropy import units
import sys

mpc_filename = sys.argv[1]
obs = EphemerisReader().read(mpc_filename)
orb = BKOrbit(obs)
print orb.summarize()
date = orb.time
for day in range(100):
    orb.predict(orb.time + 10 * units.day)
    print "{} {} {:5.2f}".format(orb.time.mpc, orb.coordinate.to_string('hmsdms', precision=2, sep=":"), orb.distance)



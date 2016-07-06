import mp_ephem
from astropy import units
import sys

mpc_filename = sys.argv[1]
obs = mp_ephem.EphemerisReader().read(mpc_filename)
orb = mp_ephem.BKOrbit(obs)
print orb.summarize()
orb.predict(orb.time + 10 * units.day)
print orb.coordinate.hmsdms, orb.distance


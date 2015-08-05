
import logging
import sys
logging.basicConfig(level=logging.INFO, format="%(filename)s %(lineno)s %(message)s")

mpc_filename = sys.argv[1]

import mpc
import orbfit_pyephem as orbfit

observations = mpc.MPCReader().read(mpc_filename)


orbit = orbfit.Orbfit(observations)

orbit.predict('2001 01 01')
print "0"*80
print "0 Position on: {} (RA DEC) -> ({})".format(orbit.date, orbit.coordinate.to_string('hmsdms', sep=':', precision=2))
print "0"*80
orbit.predict('2001 02 01')
print "0 Position on: {} (RA DEC) -> ({})".format(orbit.date, orbit.coordinate.to_string('hmsdms', sep=':', precision=2))
print "0"*80
print orbit



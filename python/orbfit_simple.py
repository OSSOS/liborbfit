__author__ = 'jjk'
"""
An example script that creates a python binding to the orbfit library.

"""
import ctypes
import tempfile
import sys


LIBORBFIT = "/usr/local/lib/liborbfit.so"


class OrbfitError(Exception):
    def __init__(self):
        super(OrbfitError, self).__init__(
            "Insufficient observations for an orbit.")


class Orbfit(object):
    """
    This class provides orbital information derived by calling 'fit_radec'.
    """

    def __init__(self, observations=None, name='', filename=None):
        """
        Given a set of observations compute the orbit using fit_radec
        and provide methods for accessing that orbit.

        Requires at least 3 observations.
        """
        self.orbfit = ctypes.CDLL(LIBORBFIT)

        if observations is None:
            observations = open(filename).readlines()
        self.observations = observations
        self.name = name
        self._fit_radec()


    @property
    def observations(self):
        return self._observations

    @observations.setter
    def observations(self, observations):
        self._observations = []
        for observation in observations:
            self._observations.append(observation.split())

    def _fit_radec(self):
        """
        call fit_radec of BK passing in the observations.

        """

        # call fitradec with mpcfile, abgfile, resfile
        self.orbfit.fitradec.restype = ctypes.POINTER(ctypes.c_double * 2)
        self.orbfit.fitradec.argtypes = [ ctypes.c_char_p, ctypes.c_char_p ]

        mpc_file = tempfile.NamedTemporaryFile(suffix='.mpc')
        for obs in self.observations:
            jd = obs[0]
            ra = obs[1]
            dec = obs[2]
            res = obs[3]
            code = obs[4]
            mpc_file.write("{} {} {} {} {}\n".format(jd, ra, dec, res, code))

        mpc_file.seek(0)

        self._abg = tempfile.NamedTemporaryFile()

        result = self.orbfit.fitradec(ctypes.c_char_p(mpc_file.name),
                                      ctypes.c_char_p(self._abg.name))

        self.distance = result.contents[0]
        self.distance_uncertainty = result.contents[1]

        self.orbfit.abg_to_aei.restype = ctypes.POINTER(ctypes.c_double * 12)
        self.orbfit.abg_to_aei.argtypes = [ ctypes.c_char_p ]
        result = self.orbfit.abg_to_aei(ctypes.c_char_p(self._abg.name))
        self.a = result.contents[0]
        self.da = result.contents[6]
        self.e = result.contents[1]
        self.de = result.contents[7]
        self.inc = result.contents[2]
        self.dinc = result.contents[8]
        self.Node = result.contents[3]
        self.dNode = result.contents[9]
        self.om = result.contents[4]
        self.dom = result.contents[10]
        self.T = result.contents[5]
        self.dT = result.contents[11]


    def __str__(self):
        """

        """
        res = "{:>10s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n".format(self.name, 
                                                            "r (AU)",
                                                            "a (AU)",
                                                            "e",
                                                            "Inc.",
                                                            "Node",
                                                            "peri.")
        res += "{:>10s} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}\n".format("fit",
                                                   self.distance,
                                                   self.a,
                                                   self.e,
                                                   self.inc,
                                                   self.Node,
                                                   self.om)
        res += "{:>10s} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}\n".format("uncert",
                                                               self.distance_uncertainty,
                                                               self.da,
                                                               self.de,
                                                               self.dinc,
                                                               self.dNode,
                                                               self.dom)

        return res

    def predict(self, jd, obs_code=568):
        """
        use the bk predict method to compute the location of the
        source on the given date.
        """
        # call predict with agbfile, jdate, obscode
        self.orbfit.predict.restype = ctypes.POINTER(ctypes.c_double * 5)
        self.orbfit.predict.argtypes = [ ctypes.c_char_p, ctypes.c_double, ctypes.c_int ]
        predict = self.orbfit.predict(ctypes.c_char_p(self._abg.name),
                                      jd,
                                      ctypes.c_int(obs_code))

        self.ra = predict.contents[0]
        self.dec = predict.contents[1]
        self.dra = predict.contents[2]
        self.ddec = predict.contents[3]
        self.pa = predict.contents[4]
        self.date = jd



if __name__ == '__main__':

    print Orbfit(filename='test_input.txt')
    

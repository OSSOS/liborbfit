from ossos import mpc

__author__ = 'jjk'

import ctypes
import tempfile
from StringIO import StringIO

import datetime
import numpy

try:
    from astropy.coordinates import ICRSCoordinates
except ImportError:
    from astropy.coordinates import ICRS as ICRSCoordinates
from astropy import units
import math

from .mpc import Time


LIBORBFIT = "/usr/local/lib/liborbfit.so"


class OrbfitError(Exception):
    def __init__(self):
        super(OrbfitError, self).__init__(
            "Insufficient observations for an orbit.")


class Orbfit(object):
    """
    This class provides orbital information derived by calling 'fit_radec'.
    """

    def __init__(self, observations):
        """
        Given a list of mpc.Observations, compute the orbit using fit_radec and provide methods for
        accessing that orbit.

        Requires at least 3 mpc.Observations in the list.
        :rtype : Orbfit
        """
        assert isinstance(observations, tuple) or isinstance(observations, list) or isinstance(observations,
                                                                                               numpy.ndarray)

        if len(observations) < 3:
            raise OrbfitError()
        self.orbfit = ctypes.CDLL(LIBORBFIT)
        self.dra = None
        self.ddec = None
        self.observations = observations
        self._fit_radec()

    def _fit_radec(self): 

        _abg_file = tempfile.NamedTemporaryFile(suffix='.abg')
        _mpc_file = tempfile.NamedTemporaryFile(suffix='.mpc')

        # call fit_radec with mpcfile and abgfile
        self.orbfit.fitradec.restype = ctypes.POINTER(ctypes.c_double * 2)
        self.orbfit.fitradec.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        for observation in self.observations:
            assert isinstance(observation, mpc.Observation)
            if observation.null_observation:
                continue
            obs = observation
            self.name = observation.provisional_name
            ra = obs.ra.replace(" ", ":")
            dec = obs.dec.replace(" ", ":")
            res = getattr(obs.comment,'plate_uncertainty',0.2)
            _mpc_file.write("{} {} {} {} {}\n".format(obs.date.jd, ra, dec, res, 568, ))
        _mpc_file.seek(0)
        result = self.orbfit.fitradec(ctypes.c_char_p(_mpc_file.name),
                                      ctypes.c_char_p(_abg_file.name))
        self.distance = result.contents[0]
        self.distance_uncertainty = result.contents[1]

        # call abg_to_aei to get elliptical elements and their chi^2 uncertainty.
        self.orbfit.abg_to_aei.restype = ctypes.POINTER(ctypes.c_double * 12)
        self.orbfit.abg_to_aei.argtypes = [ctypes.c_char_p]
        result = self.orbfit.abg_to_aei(ctypes.c_char_p(_abg_file.name))

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
        _abg_file.seek(0)
        self.abg = _abg_file.read()

    @property
    def residuals(self):
        ## compute the residuals (from the given observations)
        _residuals = ""
        for observation in self.observations:
            self.predict(observation.date)
            coord1 = ICRSCoordinates(self.coordinate.ra, self.coordinate.dec)
            coord2 = ICRSCoordinates(observation.coordinate.ra, self.coordinate.dec)
            observation.ra_residual = float(coord1.separation(coord2).arcsec)
            observation.ra_residual = (coord1.ra.degree - coord2.ra.degree) * 3600.0
            coord2 = ICRSCoordinates(self.coordinate.ra, observation.coordinate.dec)
            observation.dec_residual = float(coord1.separation(coord2).arcsec)
            observation.dec_residual = (coord1.dec.degree - coord2.dec.degree) * 3600.0
            _residuals += "{:1s}{:12s} {:+05.2f} {:+05.2f} # {}\n".format(
                observation.null_observation, observation.date, observation.ra_residual, observation.dec_residual, observation)
        return _residuals

    @property
    def arc_length(self):
        dates = []
        for observation in self.observations:
            dates.append(observation.date.jd)
        return max(dates) - min(dates)

    def __str__(self):
        """

        """

        res = "{:>10s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n".format(
            self.observations[0].provisional_name.strip(' '),
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

        res += "{:>10s} {:8.2f} days ".format("arc", self.arc_length)
        if self.dra is not None:
            res += "ephemeris uncertainty: {:8.2f} {:8.2f} {:8.2f}".format(self.dra, self.ddec, self.pa)
        res += "\n"

        return res

    def predict(self, date, obs_code=568, abg_file=None):
        """
        use the bk predict method to compute the location of the source on the given date.
        @param date: the julian date of interest or an astropy.core.time.Time object.
        @param obs_code: the Minor Planet Center observatory location code (Mauna Kea: 568 is the default)

        this methods sets the values of coordinate, dra (arcseconds), ddec (arcseconds), pa, (degrees) and date (str)
        """

        if not isinstance(date, Time):
            if isinstance(date, float):
                try:
                    date = Time(date, format='jd', scale='utc', precision=6)
                except:
                    date = None  # FIXME: this might blow up, not sure
            else:
                try:
                    date = Time(date, format='jd', scale='utc', precision=6)
                except ValueError:
                    try:
                        date = Time(date, format='mpc', scale='utc', precision=6)
                    except ValueError:
                        date = Time(date, scale='utc', precision=6)  # see if it can guess
        if hasattr(self,'time'):
	   dt = self.time - date
           if math.fabs(dt.sec) < 10:
              return
        jd = ctypes.c_double(date.jd)
        if abg_file is None:
            abg_file = tempfile.NamedTemporaryFile(suffix='.abg')
            abg_file.write(self.abg)
            abg_file.seek(0)
        # call predict with agbfile, jdate, obscode
        self.orbfit.predict.restype = ctypes.POINTER(ctypes.c_double * 5)
        self.orbfit.predict.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int]
        predict = self.orbfit.predict(ctypes.c_char_p(abg_file.name),
                                      jd,
                                      ctypes.c_int(obs_code))
        self.coordinate = ICRSCoordinates(predict.contents[0],
                                          predict.contents[1],
                                          unit=(units.degree, units.degree),
                                          obstime=date)
        self.dra = predict.contents[2]
        self.ddec = predict.contents[3]
        self.pa = predict.contents[4]
        self.date = str(date)
        self.time = date

    def rate_of_motion(self, date):
        # rate of motion at a requested date rather than averaged over the arc.
        # Date is datetime.datetime() objects.
        if isinstance(date, datetime.datetime):
            sdate = date.strftime('%Y-%m-%d')
            edate = (date + datetime.timedelta(1)).strftime('%Y-%m-%d')
        else:
            sdate = date.jd
            edate = date.jd + 1
        self.predict(edate)
        coord1 = self.coordinate
        self.predict(sdate)
        coord2 = self.coordinate
        retval = coord1.separation(coord2).arcsec / (24.)  # arcsec/hr

        return retval

    def summarize(self, date=datetime.datetime.now()):
        """Return a string summary of the orbit.

        """

        assert isinstance(date, datetime.datetime)
        at_date = date.strftime('%Y-%m-%d')
        self.predict(at_date)

        fobj = StringIO()

        # for observation in self.observations:
        # fobj.write(observation.to_string()+"\n")

        fobj.write("\n")
        fobj.write(str(self) + "\n")
        fobj.write(str(self.residuals) + "\n")
        fobj.write('arclen (days) {} '.format(self.arc_length))
        fobj.write("Expected accuracy on {:>10s}: {:6.2f}'' {:6.2f}'' moving at {:6.2f} ''/hr\n\n".format(
            at_date, self.dra, self.ddec, self.rate_of_motion(date=date)))

        fobj.seek(0)
        return fobj.read()


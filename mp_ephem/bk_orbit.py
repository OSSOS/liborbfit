import ctypes
import datetime
import logging
import math
import numpy
from StringIO import StringIO
import tempfile
import os
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.units.quantity import Quantity
from astropy.time import Time

__author__ = 'jjk'


class BKOrbitError(Exception):
    def __init__(self):
        super(BKOrbitError, self).__init__(
            "Insufficient observations for an orbit.")


class BKOrbit(object):
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
        __PATH__ = os.path.dirname(__file__)
        __ORBFIT_LIB__ = os.path.join(__PATH__, "orbit.so")
        is_64bits = ctypes.sizeof(ctypes.c_voidp) == 8

        bin_ephem = is_64bits and 'binEphem.405_64' or 'binEphem.405_32'

        os.environ['ORBIT_EPHEMERIS'] = os.getenv('ORBIT_EPHEMERIS',
                                                  os.path.join(__PATH__, 'data', bin_ephem))
        os.environ['ORBIT_OBSERVATORIES'] = os.getenv('ORBIT_OBSERVATORIES',
                                                      os.path.join(__PATH__,
                                                                   'data',
                                                                   'observatories.dat'))
        liborbfit = os.path.join(__ORBFIT_LIB__)
        self.orbfit = ctypes.CDLL(liborbfit)

        assert isinstance(observations, tuple) or isinstance(observations, list) or isinstance(observations,
                                                                                               numpy.ndarray)

        if len(observations) < 3:
            raise BKOrbitError()
        self.observations = observations
        self._a = self._e = self._inc = self._Node = self._om = self._T = None
        self._da = self._de = self._dinc = self._dNode = self._dom = self._dT = None
        self._coordinate = self._pa = self._dra = self._ddec = self._date = self._time = None
        self._fit_radec()

    def _fit_radec(self):

        _abg_file = tempfile.NamedTemporaryFile(suffix='.abg')
        _mpc_file = tempfile.NamedTemporaryFile(suffix='.mpc')

        # call fit_radec with mpcfile and abgfile
        self.orbfit.fitradec.restype = ctypes.POINTER(ctypes.c_double * 2)
        self.orbfit.fitradec.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        for observation in self.observations:
            try:
                if observation.null_observation:
                    continue
            except:
                pass
            obs = observation
            self.name = observation.provisional_name
            ra = obs.ra.replace(" ", ":")
            dec = obs.dec.replace(" ", ":")
            res = getattr(obs.comment, 'plate_uncertainty', 0.2)
            _mpc_file.write("{} {} {} {} {}\n".format(obs.date.jd, ra, dec, res, 568, ))
        _mpc_file.seek(0)
        result = self.orbfit.fitradec(ctypes.c_char_p(_mpc_file.name),
                                      ctypes.c_char_p(_abg_file.name))
        self._distance = result.contents[0] * units.AU
        self._distance_uncertainty = result.contents[1] * units.AU

        # call abg_to_aei to get elliptical elements and their chi^2 uncertainty.
        self.orbfit.abg_to_aei.restype = ctypes.POINTER(ctypes.c_double * 13)
        self.orbfit.abg_to_aei.argtypes = [ctypes.c_char_p]
        result = self.orbfit.abg_to_aei(ctypes.c_char_p(_abg_file.name))

        self._a = result.contents[0] * units.AU
        self._da = result.contents[6] * units.AU
        self._e = result.contents[1] * units.dimensionless_unscaled
        self._de = result.contents[7] * units.dimensionless_unscaled
        self._inc = result.contents[2] * units.degree
        self._dinc = result.contents[8] * units.degree
        self._Node = result.contents[3] * units.degree
        self._dNode = result.contents[9] * units.degree
        self._om = result.contents[4] * units.degree
        self._dom = result.contents[10] * units.degree
        self._T = result.contents[5] * units.day
        self._dT = result.contents[11] * units.day
        self._epoch = Time(result.contents[12] * units.day, scale='utc', format='jd')
        _abg_file.seek(0)
        self.abg = _abg_file.read()

    @property
    def a(self):
        """
        :return: semi-major axis.
        :rtype: Quantity
        """
        return self._a

    @property
    def epoch(self):
        """
        :return: the epoch of the orbital elements.
        :rtype: Quantity
        """
        return self._epoch

    @property
    def da(self):
        """
        :return:Uncertainty of orbital semi-major axis.
        :rtype: Quantity
        """
        return self._da

    @property
    def e(self):
        """
        :return: eccentricity.
        :rtype: Quantity
        """
        return self._e

    @property
    def de(self):
        """
        :return: of orbital eccentricity
        :rtype: Quantity
        """
        return self._de

    @property
    def inc(self):
        """
        :return: Inclination
        :rtype: Quantity
        """
        return self._inc

    @property
    def dinc(self):
        """

        :return: Uncertainty of orbital Inclination
        :rtype: Quantity
        """
        return self._dinc

    @property
    def om(self):
        """

        :return: Argument of peri-centre
        :rtype: Quantity
        """
        return self._om

    @property
    def dom(self):
        """
        :rtype : Quantity
        :return: Uncertainty in argument of peri-centre
        """
        return self._dom

    @property
    def Node(self):
        """
        :rtype : Quantity
        :return: Ascending Node
        """
        return self._Node

    @property
    def dNode(self):
        """
        :rtype : Quantity
        :return: The uncertainty in the chi2 fit of ascending Node value.
        """
        return self._dNode

    @property
    def T(self):
        """
        :return: Time of pericentre passage
        :rtype: Quantity
        """
        return self._T

    @property
    def dT(self):
        """
        :return: Uncertainty in time of pericentre passage
        :rtype: Quantity
        """
        return self._dT

    @property
    def residuals(self, overall=False):
        """
        Builds a summary of the residuals of a fit.  This is useful for visually examining the
        goodness of fit for a small number of observations and the impact of adding a few tentative observations.
        :param overall: should we return the string with the residuals or a list containing the overall residuals?
        :return: A string representation of the residuals between the best fit orbit and  astrometric measurements.
        :rtype: str
        """
        _residuals = ""
        overall_resids = []
        for observation in self.observations:
            self.predict(observation.date)
            coord1 = SkyCoord(self.coordinate.ra, self.coordinate.dec)
            coord2 = SkyCoord(observation.coordinate.ra, self.coordinate.dec)
            observation.ra_residual = float(coord1.separation(coord2).arcsec)
            observation.ra_residual = (coord1.ra.degree - coord2.ra.degree) * 3600.0
            coord2 = SkyCoord(self.coordinate.ra, observation.coordinate.dec)
            observation.dec_residual = float(coord1.separation(coord2).arcsec)
            observation.dec_residual = (coord1.dec.degree - coord2.dec.degree) * 3600.0
            overall_resids.append(math.sqrt(observation.ra_residual ** 2 + observation.dec_residual ** 2))
            _residuals += "{:1s}{:12s} {:+05.2f} {:+05.2f} # {}\n".format(
                observation.null_observation, observation.date, observation.ra_residual, observation.dec_residual,
                observation)
        if overall:
            return overall_resids  # values in arcsec
        else:
            return _residuals

    @property
    def coordinate(self):
        """
        A coordinate object that holds the current predicted location of the object. (see self.time/self.date)
        :return: location of the source determined by predict
        :rtype: ICRS
        """
        return self._coordinate

    @property
    def dra(self):
        """
        :return: Major axis of predicted location uncertainty ellipse
        :rtype: Quantity
        """
        return self._dra

    @property
    def ddec(self):
        """
        :return: Minor axis of predicted location uncertainty ellipse
        :rtype: Quantity
        """
        return self._ddec

    @property
    def date(self):
        """
        :return: String representation of the date prediction is for.
        :rtype: str
        """
        return self._date

    @property
    def time(self):
        """
        :return: Time the predicted location was computed for.
        :rtype: time.Time
        """
        return self._time

    @property
    def pa(self):
        """
        :return: The position angle of the uncertainty of the prediction ellipse.
        :rtype: Quantity
        """
        return self._pa

    @property
    def arc_length(self):
        """
        :return: the arc length of the observations used in the current fit of the orbit.
        :rtype: Quantity
        """
        dates = []
        for observation in self.observations:
            dates.append(observation.date.jd)
        return (max(dates) - min(dates)) * units.day

    @property
    def distance(self):
        """
        Estimate of the current geocentric distance to the source.
        :return:
        """
        return self._distance

    @property
    def distance_uncertainty(self):
        """
        Estimate of the uncertainty in the geocentric distance estimate.
        :return:
        """
        return self._distance_uncertainty

    def __str__(self):
        """

        """

        res = "{:>12s} {:>10s} {:>10s} {:>7s} {:>10s} {:>13s} {:>10s} on JD {:15}\n".format(
            self.observations[0].provisional_name.strip(' '),
            "distance",
            "a ",
            "e",
            "Inc.",
            "Node",
            "peri.",
            self.epoch.iso)
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

        res += "{:>10s} {:8.2f} ".format("arc", self.arc_length)
        if self.dra is not None:
            res += "ephemeris uncertainty: {:8.2f} {:8.2f} {:8.2f}".format(self.dra, self.ddec, self.pa)
        res += "\n"

        return res

    def predict(self, date, obs_code=568, abg_file=None, minimum_delta=None):
        """
        use the bk predict method to compute the location of the source on the given date.
        @param date: the julian date of interest or an astropy.core.time.Time object.
        @param obs_code: the Minor Planet Center observatory location code (Mauna Kea: 568 is the default)

        this methods sets the values of coordinate, dra (arcseconds), ddec (arcseconds), pa, (degrees) and date (str)
        :param abg_file: name of the file containing the BK orbit ABG data.
        :param minimum_delta: the smallest time offset between computing positions, otherwise use previous position.

        """

        if minimum_delta is None:
            minimum_delta = 0.00001 * units.day

        if not isinstance(date, Time):
            if isinstance(date, float):
                try:
                    date = Time(date, format='jd', scale='utc', precision=6)
                except:
                    raise ValueError("Bad date value: {}".format(date))
            date = Time(date)

        # for speed reasons we only compute positions every 10 seconds.
        if hasattr(self, 'time') and isinstance(self.time, Time):
            if -minimum_delta < self.time - date < minimum_delta:
                return

        jd = ctypes.c_double(date.jd)
        if abg_file is None:
            abg_file = tempfile.NamedTemporaryFile(suffix='.abg')
            abg_file.write(self.abg)
            abg_file.seek(0)

        # call predict with agbfile, jdate, obscode
        self.orbfit.predict.restype = ctypes.POINTER(ctypes.c_double * 6)
        self.orbfit.predict.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int]
        predict = self.orbfit.predict(ctypes.c_char_p(abg_file.name),
                                      jd,
                                      ctypes.c_int(obs_code))
        self._coordinate = SkyCoord(predict.contents[0],
                                    predict.contents[1],
                                    unit=(units.degree, units.degree),
                                    obstime=date)
        self._dra = predict.contents[2] * units.arcsec
        self._ddec = predict.contents[3] * units.arcsec
        self._pa = predict.contents[4] * units.degree
        self._distance = predict.contents[5] * units.AU
        self._date = str(date)
        self._time = date

    def rate_of_motion(self, date=None):
        """
        The rate of motion of the date of the current predicted location, or the date given if specified.
        :param date: Time
        :return: rate of motion
        :rtype: Quantity
        """
        # rate of motion at a requested date rather than averaged over the arc.
        # Date is datetime.datetime() objects.
        if date is None:
            if self.time is None:
                date = datetime.datetime.now()
            else:
                date = self.time
        if date is None:
            return None
        if isinstance(date, datetime.datetime):
            try:
                date = date.strftime('%Y-%m-%d %H:%M:%S')
            except Exception as e:
                logging.error(str(e))
                return None
        if isinstance(date, str):
            date = Time(date, scale='utc')

        try:
            sdate = date.jd
            edate = date.jd + 1
        except Exception as e:
            logging.error(str(e))
            return None
        self.predict(edate)
        coord1 = self.coordinate
        self.predict(sdate)
        coord2 = self.coordinate
        return coord1.separation(coord2).to(units.arcsec) / (24. * units.hour)  # arcsec/hr

    def summarize(self, date=None):
        """Return a string summary of the orbit.
        :param date: Date to make the summary for
        :rtype: str
        """

        if date is None:
            date = datetime.datetime.now()

        at_date = Time(date)
        self.predict(at_date)

        fobj = StringIO()

        # for observation in self.observations:
        # fobj.write(observation.to_string()+"\n")

        fobj.write("\n")
        fobj.write(str(self) + "\n")
        fobj.write(str(self.residuals) + "\n")
        fobj.write('arclen {} '.format(self.arc_length))
        fobj.write("Expected accuracy on {:>10s}: {:6.2f} {:6.2f} moving at {:6.2f} \n\n".format(
            at_date, self.dra, self.ddec, self.rate_of_motion(date=date)))

        fobj.seek(0)
        return fobj.read()
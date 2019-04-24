import ctypes
import datetime
import glob
import logging
import math
import numpy
import tempfile
import os
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.units.quantity import Quantity
from astropy.time import Time
from .ephem import EphemerisReader, ObsRecord

__author__ = 'jjk'


class BKOrbitError(Exception):
    def __init__(self):
        super(BKOrbitError, self).__init__(
            "Insufficient observations for an orbit.")


class BKOrbit(object):
    """
    This class provides orbital information derived by calling 'fit_radec'.
    """

    def __init__(self, observations, ast_filename=None, abg_file=None):
        """
        Given a list of mpc.Observations, compute the orbit using fit_radec and provide methods for
        accessing that orbit.

        Requires at least 3 mpc.Observations in the list.
        :rtype : Orbfit
        """
        __PATH__ = os.path.dirname(__file__)
        # find the orbit.so library
        __ORBFIT_LIB__ = glob.glob(os.path.join(__PATH__, '*orbit*.so'))[0]

        # Choose a format of JPL binary ephemeris file
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

        if ast_filename is None:
            assert isinstance(observations, tuple) or isinstance(observations, list) or isinstance(observations,
                                                                                                   numpy.ndarray)

        self._observations = observations
        self._r_mag = None
        self.abg_filename = abg_file
        self.ast_filename = ast_filename
        self.abg = None
        self._a = self._e = self._inc = self._Node = self._om = self._T = None
        self._da = self._de = self._dinc = self._dNode = self._dom = self._dT = None
        self._coordinate = self._pa = self._dra = self._ddec = self._date = self._time = None
        self._distance = None
        self._distance_uncertainty = None
        self._residuals = None
        self._overall_residuals = []
        self._fit_radec()

    @property
    def observations(self):
        if self._observations is None:
            self._observations = EphemerisReader().read(self.ast_filename)
        return self._observations

    def compute_median_mag(self, band):
        mags = []
        for observation in self.observations:
            if observation.mag is not None and observation.band is not None and observation.band.lower() == band:
                mags.append(float(observation.mag))
        if len(mags) > 0:
            return numpy.percentile(mags, 50)
        return None

    @property
    def r_mag(self):
        offsets = {'r': 0.0, 'g': -0.7, 'i': +0.4, 'w': -0.2}
        bands = ['r', 'g', 'i', 'w']
        if self._r_mag is None:
            for band in bands:
                mag = self.compute_median_mag(band)
                if mag is not None:
                    self._r_mag = mag + offsets[band]
                    break
        return self._r_mag
    
    @property
    def mag(self):
        return self.r_mag

    def _fit_radec(self):
        # call fit_radec with mpcfile and abgfile
        self.orbfit.fitradec.restype = ctypes.POINTER(ctypes.c_double * 2)
        self.orbfit.fitradec.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        build_abg = True
        if self.abg_filename is not None and os.access(self.abg_filename, os.R_OK):
            build_abg = False
            _abg_file = open(self.abg_filename, 'r')
        if build_abg:
            _mpc_file = tempfile.NamedTemporaryFile(mode='w', suffix='.mpc')
            if len(self.observations) < 3:
                raise BKOrbitError()

            for observation in self.observations:
                try:
                    if observation.null_observation:
                        continue
                except Exception:
                    pass
                assert isinstance(observation, ObsRecord)
                obs = observation
                self.name = obs.provisional_name
                ra = obs.ra.replace(" ", ":")
                dec = obs.dec.replace(" ", ":")
                res = getattr(obs.comment, 'plate_uncertainty', 0.2)
                if obs.location is None:
                    _mpc_file.write("{} {} {} {} {}\n".format(obs.date.jd, ra, dec, res, int(obs.observatory_code) ))
                else:
                    _mpc_file.write("{} {} {} {} {} {} {}\n".format(obs.date.jd, ra, dec, res,
                                                                    obs.location.x,
                                                                    obs.location.y,
                                                                    obs.location.z))

            _mpc_file.seek(0)
            if self.abg_filename is None:
                _abg_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.abg')
                _abg_file_name = _abg_file.name
            else:
                _abg_file = open(self.abg_filename, 'w+')
                _abg_file_name = _abg_file.name
            result = self.orbfit.fitradec(ctypes.c_char_p(bytes(_mpc_file.name, 'utf-8')),
                                          ctypes.c_char_p(bytes(_abg_file_name, 'utf-8')))

        _abg_file.seek(0)

        # call abg_to_aei to get elliptical elements and their chi^2 uncertainty.
        self.orbfit.abg_to_aei.restype = ctypes.POINTER(ctypes.c_double * 15)
        self.orbfit.abg_to_aei.argtypes = [ctypes.c_char_p]
        result = self.orbfit.abg_to_aei(ctypes.c_char_p(bytes(_abg_file.name, 'utf-8')))
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
        self._distance = result.contents[13] * units.au
        self._distance_uncertainty = result.contents[14] * units.au
        self.abg = _abg_file.read()
        self._residuals = None

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
    def ra(self):
        """
        :return: Right Ascension in Equatorial J2000 coordinates.
        :rtype: Quantity
        """
        return self.coordinate.ra

    @property
    def dec(self):
        """
        :return: Declination in Equatorial J2000 coordinates.
        :rtype: Quantity
        """
        return self.coordinate.dec

    def compute_residuals(self):
        """
        Builds a summary of the residuals of a fit and loads those into the observation objects.
        """
        for observation in self.observations:
            self.predict(observation.date, obs_code=int(observation.observatory_code))
            coord1 = SkyCoord(self.coordinate.ra, self.coordinate.dec)
            coord2 = SkyCoord(observation.coordinate.ra, self.coordinate.dec)
            # observation.ra_residual = float(coord1.separation(coord2).arcsec)
            observation.ra_residual = (coord1.ra.degree - coord2.ra.degree) * 3600.0
            coord2 = SkyCoord(self.coordinate.ra, observation.coordinate.dec)
            # observation.dec_residual = float(coord1.separation(coord2).arcsec)
            observation.dec_residual = (coord1.dec.degree - coord2.dec.degree) * 3600.0

    @property
    def overall_residuals(self):
        _overall_residuals = []
        for observation in self.observations:
            _overall_residuals.append(math.sqrt(observation.ra_residual ** 2 + observation.dec_residual ** 2))
        return _overall_residuals  # values in arcsec

    @property
    def residuals(self):
        if self._residuals is None:
            self.compute_residuals()
            _residuals = ""
            for observation in self.observations:
                _residuals += "{:1s}{:12s} {:+05.2f} {:+05.2f} # {}\n".format(
                    str(observation.null_observation),
                    str(observation.date), observation.ra_residual,
                    observation.dec_residual,
                    observation)
            self._residuals = _residuals
        return self._residuals

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

    def __format__(self, format_spec):
        return self.__str__()

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
        res += "{:>10s} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}\n".format("uncertainty",
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

        this methods sets the values of coordinate, dra (arc seconds), ddec (arc seconds), pa, (degrees) and date (str)

        @param date: the julian date of interest or an astropy.core.time.Time object.
        @param obs_code: the Minor Planet Center observatory location code (Maunakea: 568 is the default)
        @param minimum_delta: minimum difference in time between recomputing a predicted location.
        @param abg_file: the 'abg' formatted file to use for the prediction.
        """

        if minimum_delta is None:
            # Set the minimum delta to very small value so that we always compute new positions
            minimum_delta = 0.00001 * units.second

        if not isinstance(date, Time):
            if isinstance(date, float):
                try:
                    date = Time(date, format='jd', scale='utc', precision=6)
                except:
                    raise ValueError("Bad date value: {}".format(date))
            _date = Time(date)
        else:
            _date = date

        # for speed reasons we only compute positions every 10 seconds.
        if hasattr(self, 'time') and isinstance(self.time, Time):
            if -minimum_delta < self.time - _date < minimum_delta:
                return

        jd = ctypes.c_double(_date.jd)
        if abg_file is None:
            abg_file = tempfile.NamedTemporaryFile(suffix='.abg')
            abg_file.write(bytes(self.abg, 'utf-8'))
            abg_file.seek(0)

        self.orbfit.predict.restype = ctypes.POINTER(ctypes.c_double * 6)
        self.orbfit.predict.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int]
        predict = self.orbfit.predict(ctypes.c_char_p(bytes(abg_file.name, 'utf-8')),
                                      jd,
                                      ctypes.c_int(obs_code))
        self._coordinate = SkyCoord(predict.contents[0],
                                    predict.contents[1],
                                    unit=(units.degree, units.degree),
                                    obstime=_date)
        self._dra = predict.contents[2] * units.arcsec
        self._ddec = predict.contents[3] * units.arcsec
        self._pa = predict.contents[4] * units.degree
        self._distance = predict.contents[5] * units.AU
        self._date = str(_date)
        self._time = _date

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
            start_date = date.jd
            end_date = date.jd + 1
        except Exception as e:
            logging.error(str(e))
            return None
        self.predict(end_date)
        coord1 = self.coordinate
        self.predict(start_date)
        coord2 = self.coordinate
        return coord1.separation(coord2).to(units.arcsec) / (24. * units.hour)  # arc seconds / hr

    def summarize(self, date=None):
        """Return a string summary of the orbit.
        :param date: Date to make the summary for
        :rtype: str
        """

        if date is None:
            date = datetime.datetime.now()

        at_date = Time(date)
        self.predict(at_date)

        summary = "\n"
        summary += str(self) + "\n"
        summary += str(self.residuals) + "\n"
        summary += 'arc length {} '.format(self.arc_length)
        summary += "Expected accuracy on {:>10s}: {:6.2f} {:6.2f} moving at {:6.2f} \n\n".format(
            str(at_date), self.dra, self.ddec, self.rate_of_motion(date=date))

        return summary

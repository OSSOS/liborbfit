"""
send queries to the Horizons batch query system via a web UI.  Derived from code supplied by
Wes Fraser and Michele Bannister with more smarty pants stuff added by JJ Kavelaars.
            irint self.datadd
"""
import copy
import re
import requests
import scipy
import logging

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units
from astropy.units.quantity import Quantity
from astropy.io.ascii import Csv
from astropy.table import Table


class Query(object):
    """
    A query object is used to build the query structure that will be sent to Horizons web system.
    """

    SERVER = 'ssd.jpl.nasa.gov'
    PROTOCOL = 'https'
    END_POINT = 'horizons_batch.cgi'

    horizons_quantities = [0, 'Astrometric RA & DEC', 'Apparent RA & DEC',
                           'Rates; RA & DEC', 'Apparent AZ & EL', 'Rates; AZ & EL',
                           'Sat. X & Y, pos. ang', 'Local app. sid. time',
                           'Airmass', 'APmag', 'Illuminated fraction',
                           'Defect of illumin.', 'Sat. angle separ/vis',
                           'Target angular diam.', 'Obs sub-lng & sub-lat',
                           'Sun sub-long & sub-lat', 'Sub Sun Pos. Ang & Dis',
                           'N. Pole Pos. Ang & Dis', 'Helio eclip. lon & lat',
                           'Helio range & rng rate', 'Obsrv range & rng rate',
                           'One-Way Light-Time', 'Speed wrt Sun & obsrvr',
                           'Sun-Obs-Targ ELONG ang', 'Sun-Targ-Obs PHASE ang',
                           'Targ-Obsrv-Moon/Illum%', 'Obs-Primary-Targ angl',
                           'Pos. Ang;radius & -vel', 'Orbit plane angle',
                           'Constellation ID', 'Delta-T (CT - UT)', 'Obs eclip. lon & lat',
                           'North pole RA & DEC', 'Galactic latitude', 'Local app. SOLAR time',
                           'Earth->Site lt-time', 'RA & DEC uncertainty', 'POS error ellipse',
                           'POS uncertainty (RSS)', 'Range & Rng-rate sig.', 'Doppler/delay sigmas',
                           'True anomaly angle', 'Local app. hour angle', 'PHASE angle & bisector']

    default_quantities = ['Astrometric RA & DEC', 'Rates; RA & DEC',
                          'RA & DEC uncertainty', 'APmag',
                          'POS error ellipse', 'Helio range & rng rate', 'Obsrv range & rng rate',
                          'PHASE angle & bisector', 'Sun-Obs-Targ ELONG ang']

    default_horizons_params = {'batch': "'{}'".format(1),
                               'COMMAND': "'{}'".format('Ceres'),
                               'MAKE_EPHEM': "'YES'",
                               'TABLE_TYPE': "'OBSERVER'",
                               'CAL_FORMAT': "'BOTH'",
                               'TIME_DIGITS': "'SECONDS'",
                               'ANG_FORMAT': "'DEG'",
                               'CENTER': "'568@399'",
                               'START_TIME': "'JD {}'".format(Time.now()),
                               'STOP_TIME': "'JD {}'".format(Time.now() + 10*units.day),
                               'STEP_SIZE': "'{} {}'".format(1, 'd'),
                               'QUANTITIES': None,
                               'REF_SYSTEM': "'J2000'",
                               'SKIP_DAYLT': "'NO'",
                               'EXTRA_PREC': "'NO'",
                               'R_T_S_ONLY': "'NO'",
                               'CSV_FORMAT': "'YES'"}

    def __init__(self, target, start_time, stop_time, step_size):
        """

        @type start_time: Time
        @type stop_time: Time
        @type step_size: Quantity
        @type target: str
        @return: Query
        """
        self._target = target
        self._params = copy.copy(Query.default_horizons_params)
        self._quantities = None
        self.quantities = Query.default_quantities
        self._data = None
        self._start_time = start_time
        self._stop_time = stop_time
        self._step_size = step_size
        self._center = "568"
        self._nobs = None
        self._arc_length = None
        self._current_time = None

    @property
    def target(self):
        """
        The solar system body to get an ephemeris for.

        @return: str
        """
        return self._target

    @target.setter
    def target(self, target):
        if self._target != target:
            self._target = target
            self.reset()

    def reset(self):
        self._data = None

    @property
    def center(self):
        """
        The location of the observer.

        @return: str
        """
        return self._center

    @center.setter
    def center(self, center):
        if "@" in str(center):
            self._center = center
        else:
            self._center = "{}@399".format(center)

    @property
    def command(self):
        """
        The command to send to horizons, this is the same as the targeted but formatted for Horizons
        @return: str
        """
        return "'{}'".format(self.target)

    @property
    def params(self):
        """
        Build the params set for the Horizons query.  There is a default set that is overridden by matching attributes.
        @return: dict
        """
        for key in self._params:
            attr = key.lower()
            if hasattr(self, attr):
                self._params[key] = getattr(self, attr, self._params[key])
        return self._params

    @params.setter
    def params(self, **kwargs):
        """
        Set the parameters to be set on the Horizons query.

        See Query.default_horizons_params for a list of possible parameters.
        @return:
        """
        for (key, value) in enumerate(kwargs):
            self._params[key] = value

        # this is a required setting, or we have trouble with parsing.
        self._params['CSV_FORMAT'] = "'YES'"

    @property
    def quantities(self):
        """
        An array of quantities to request from Horizons for the ephemeris.  Returned as a string for Horizons query.

        @return: str
        """
        s = "{}".format(self._quantities)
        return "'{}'".format(s[1:-1])

    @quantities.setter
    def quantities(self, quantities=None):
        if quantities is None:
            quantities = []
        quantities.extend(Query.default_quantities)

        self._quantities = []
        for quantity in quantities:
            try:
                idx = int(quantity)
            except:
                idx = Query.horizons_quantities.index(quantity)
            if idx is not None and idx not in self._quantities:
                self._quantities.append(idx)

    @property
    def start_time(self):
        """
        Start time of the ephemeris, expressed as Julian Date and in string formatted expected by Horizons.

        This attribute is set to a Time object but returned as a Horizons formatted string, same for stop_time.

        @return: str
        """
        return "'JD {}'".format(self._start_time.jd)

    @start_time.setter
    def start_time(self, start_time):
        assert isinstance(start_time, Time)
        if start_time != self._start_time:
            self.reset()
            self._start_time = start_time

    @property
    def stop_time(self):
        """
        Stop time of the ephemeris, expressed as Julian Date and in string formatted expected by Horizons.

        This attribute is set to a Time object but returned as a Horizons formatted string, same for start_time.

        @return: str
        """
        return "'JD {}'".format(self._stop_time.jd)

    @stop_time.setter
    def stop_time(self, stop_time):
        assert isinstance(stop_time, Time)
        if self._stop_time != stop_time:
            self._stop_time = stop_time
            self.reset()

    @property
    def step_size(self):
        """
        Size of the step in the ephemeris.

        This is set as a Quantity object with time dimension but returned as a string in Horizons query format.

        @return: str.
        """
        s = "{:1.0f}".format(self._step_size)[:3]
        return "'{}'".format(s)

    @step_size.setter
    def step_size(self, step_size):
        assert isinstance(step_size, Quantity)
        if self._step_size != step_size:
            self._step_size = step_size
            self.reset()

    @property
    def data(self):
        """
        The a list of strings returned from the Horizons query.

        @return: list
        """
        if self._data is None:
            self.query()
        return self._data

    def query(self):
        """
        Connect to the service and make the query.

        @return:
        """
        url = '{}://{}/{}'.format(Query.PROTOCOL,
                                  Query.SERVER,
                                  Query.END_POINT)
        logging.info("Sending JPL/Hoirzons query.\n")
        response = requests.get(url, params=self.params)
        response.raise_for_status()
        self._data = []
        for line in response.iter_lines():
            self._data.append(line)


class Body(object):
    """An Horizons Ephemeris returned as a result of a query to the JPL/Horizons.

    """

    def __init__(self, name, start_time=None, stop_time=None, step_size=None, center=None):
        """

        @rtype: Ephemeris
        @param name: Body to build an ephemeris for
        @param start_time:  start time of the ephemeris
        @type start_time: Time
        @param stop_time: stop time of the ephemeris
        @type stop_time: Time
        @param step_size: size of time step for ephemeris
        @type step_size: Quantity
        """
        self.name = str(name)

        # make sure the input quantities are reasonable.
        if start_time is None:
            start_time = Time.now()
        self._start_time = Time(start_time, scale='utc')

        if stop_time is None:
            stop_time = Time.now() + 1.0 * units.day
        self._stop_time = Time(stop_time, scale='utc')
        step_size = step_size is None and 1 or step_size

        if not isinstance(step_size, Quantity):
            step_size *= units.day
        self.step_size = step_size
        if center is None:
            center = 568
        self._center = center
        self._ephemeris = None
        self._elements = None
        self._current_time = None
        self._data = None

    def __str__(self):
        s = '{:>5}: {}\n'.format('Name', self.name)
        for attr in ['a', 'e', 'Inc', 'Omega', 'omega', 'M', 'Epoch', 'n']:
            s += "{:>5}: {:<12}\n".format(attr, getattr(self, attr))
        return s

    def _reset(self):
        self._ephemeris = None
        self._elements = None
        self._nobs = None
        self._arc_length = None
        self._data = None

    def __call__(self, keyword, **kwargs):
        if keyword not in self.ephemeris.colnames:
            raise KeyError("Requested column {} does not in: {}".format(keyword, self.ephemeris.colnames))
        return scipy.interp(self.current_time.jd,
                            self.ephemeris['Time'].jd,
                            self.ephemeris[keyword])

    @property
    def elongation(self):
        """
        Sun-Observer-Target Elongation
        S-O-T /r
        :return:
        """
        return scipy.interp(self.current_time.jd,
                            self.ephemeris['Time'].jd,
                            self.ephemeris['S-O-T']) * units.degree


    @property
    def start_time(self):
        """
        Start time of the ephemeris that was retrieved from Horizons.
        @rtype:  Time
        """
        if self._start_time is None:
            self._start_time = min(self.ephemeris['Time'])
        return self._start_time

    @property
    def stop_time(self):
        """
        Stop time of the ephemeris that was retrieved from Horizons.
        @rtype: Time
        """
        if self._stop_time is None:
            self._stop_time = max(self.ephemeris['Time'])
        return self._stop_time

    def _parse_ephemeris(self):
        """Parse the ephemeris out of the responses from Horizons and place in self._ephemeris."""

        start_of_ephemeris = b'$$SOE'
        end_of_ephemeris = b'$$EOE'
        start_of_failure = b'!$$SOF'
        start_idx = None
        end_idx = None
        for idx, line in enumerate(self.data):
            if line.startswith(start_of_ephemeris):
                start_idx = idx
            if line.startswith(end_of_ephemeris):
                end_idx = idx

        if start_idx is None or end_idx is None:
            fail_idx = None
            for idx, line in enumerate(self.data):
                if start_of_failure in line:
                    fail_idx = idx
                    break
            msg = fail_idx is None and self.data or self.data[fail_idx-2]
            logging.error(msg)
            raise ValueError(msg, "failed to build ephemeris")

        # the header of the CSV structure is 2 lines before the start_of_ephmeris
        csv_lines = [str(self.data[start_idx - 2], 'utf-8')]
        for l in self.data[start_idx+1:end_idx]:
            line = str(l, 'utf-8')
            csv_lines.append(line)
        csv = Csv()
        table = csv.read(csv_lines)
        try:
            table['Time'] = Time(table['Date_________JDUT'], format='jd')
        except KeyError:
            raise ValueError(self.data, "Horizons result did not contain a JD Time column, rebuild the Query.")
        self._ephemeris = table

    def _parse_elements(self):
        """Parse the elements out of the response from Horizons and place in elements object."""

        # the elements are in a '****' record at the start of file
        elements_record_started = False
        self._elements = {}
        for line in self.data:
            if line.startswith('****'):
                if elements_record_started:
                    break
                elements_record_started = True
                continue
            for part in re.findall('((\S+=)\s+(\S+))', line):
                key = part[1].strip().strip('=')
                try:
                    value = float(part[2].strip())
                except:
                    value = part[2].strip()
                self._elements[key] = value

    def _parse_obs_arc(self):
        parts = re.search('# obs: (\d+) \((\d+)-(\d+)\)', str(self.data))
        if parts is None:
            # Just fake some data
            self._arc_length = 30 * units.day
            self._nobs = 3
        else:
            self._arc_length = (int(parts.group(3)) - int(parts.group(2))) * units.year
            self._nobs = int(parts.group(1))

    @property
    def data(self):
        """

        :rtype: list
        """
        if self._data is None:
            q = Query(self.name, self._start_time, self._stop_time, self.step_size)
            q.center = self._center
            self._data = q.data
        return self._data

    @property
    def elements(self):
        """
        The orbital elements for the body, as returned by JPL/Horizons.

        @rtype: dict
        """
        if self._elements is None:
            self._parse_elements()
        return self._elements

    @property
    def ephemeris(self):
        """
        A table containing the ephemeris of the body.

        @rtype: Table
        """
        if self._ephemeris is None:
            self._parse_ephemeris()
        return self._ephemeris

    @property
    def nobs(self):
        """
        The number of observations used to determine the orbit used to build the ephmeris.

        @rtype: int
        """
        if self._nobs is None:
            self._parse_obs_arc()
        return self._nobs

    @property
    def arc_length(self):
        """
        The length of the observed arc used to determine the orbit used to build the ephemeris.

        @rtype: Quantity
        """
        if self._arc_length is None:
            self._parse_obs_arc()
        return self._arc_length

    @property
    def current_time(self):
        """
        Current time for position predictions.

        If current time is set to a value outside the bounds of the available ephemeris a new query to JPL/Horizons
        occurs.

        @rtype: Time
        @return: the time of the current ra/dec/rates selected from the ephmeris.
        """
        if self._current_time is None:
            self.current_time = self.stop_time + (self.stop_time - self.start_time)/2.0
        return self._current_time

    @current_time.setter
    def current_time(self, current_time):
        self._current_time = Time(current_time, scale='utc')
        if not (self.stop_time >= self.current_time >= self.start_time):
            logging.info("Resetting the ephemeris time boundaries")
            self._start_time = Time(self.current_time - 10.0*units.minute)
            self._stop_time = Time(self.current_time + 10.0*units.minute)
            self.step_size = 5*units.minute
            self._reset()

    @property
    def coordinate(self):
        """
        Prediction position of the target at current_time

        @rtype: SkyCoord
        """

        ra = scipy.interp(self.current_time.jd,
                          self.ephemeris['Time'].jd,
                          self.ephemeris['R.A._(ICRF/J2000.0)']) * units.degree
        dec = scipy.interp(self.current_time.jd,
                           self.ephemeris['Time'].jd,
                           self.ephemeris['DEC_(ICRF/J2000.0)']) * units.degree
        distance = 40 * units.au
        return SkyCoord(ra, dec, distance=distance)

    @property
    def ra_rate(self):
        """
        Uncertainty in the prediction location.

        @rtype: Quantity angle/time
        """

        ra_rate = scipy.interp(self.current_time.jd,
                               self.ephemeris['Time'].jd,
                               self.ephemeris['dRA*cosD']) * units.arcsec / units.hour

        return ra_rate

    @property
    def dec_rate(self):
        """
        Uncertainty in the prediction location.

        @rtpye: Quantity  angle/time
        """

        dec_rate = scipy.interp(self.current_time.jd,
                                self.ephemeris['Time'].jd,
                                self.ephemeris['d(DEC)/dt']) * units.arcsec / units.hour

        return dec_rate

    @property
    def dra(self):
        """
        Uncertainty in the prediction location.

        @rtpye: Quantity  angle
        """

        dra = scipy.interp(self.current_time.jd,
                           self.ephemeris['Time'].jd,
                           self.ephemeris['RA_3sigma']) * units.arcsec

        return dra

    @property
    def mag(self):
        """
        Visual magnitude of the source.
        """
        return scipy.interp(self.current_time.jd,
                            self.ephemeris['Time'].jd,
                            self.ephemeris['APmag'])

    @property
    def alpha(self):
        """
        Phase angle.
        """
        return scipy.interp(self.current_time.jd,
                            self.ephemeris['Time'].jd,
                            self.ephemeris['phi']) * units.degree

    @property
    def ddec(self):
        """
        Uncertainty in the prediction location.

        @rtype: Quantity  angle
        """
        return scipy.interp(self.current_time.jd,
                            self.ephemeris['Time'].jd,
                            self.ephemeris['DEC_3sigma']) * units.arcsec

    @property
    def pa(self):
        """
        Plane Of Sky angle of the ra/dec uncertainty ellipse.

        @rtype: Quantity
        """
        return scipy.interp(self.current_time.jd,
                            self.ephemeris['Time'].jd,
                            self.ephemeris['Theta']) * units.degree

    def predict(self, current_time):
        """
        The time for which coordinate and other attributes will be computed.

        This function exists so that a horizons.Ephemeris object can be used where a orbfit.Orbfit object is used.
        @type current_time: Time
        @param current_time: Time
        """
        self.current_time = current_time

    @property
    def a(self):
        """
        Orbital semi-major axes
        @return: Quantity AU

        """
        return self.elements['A'] * units.AU

    @property
    def e(self):
        """
        orbital eccentricity

        @return: float
        """
        return self.elements['EC']

    @property
    def Inc(self):
        """
        orbital inclination

        @return: Quantity degree
        """
        return self.elements['IN'] * units.degree

    @property
    def M(self):
        """
        orbital Mean anomaly

        @return: Quantity degree
        """
        return self.elements['MA'] * units.degree

    @property
    def Omega(self):
        """
        orbital argument of ascending node.

        @return: Quantity degree
        """
        return self.elements['OM'] * units.degree

    @property
    def omega(self):
        """
        orbital argument of peri-focus
        @return:
        """
        return self.elements['W'] * units.degree

    @property
    def Epoch(self):
        """
        epoch of the orbital elements  in dynamical barycenter units
        @return: Time
        """
        return Time(self.elements.get('EPOCH', Time.now().jd), format='jd')

    @property
    def n(self):
        """
        mean motion
        @return: Quantity degree / day
        """
        return self.elements['N'] * units.degree / units.day

    @property
    def delta(self):
        """
        Distance from Observer to Target
        :return:
        """
        return scipy.interp(self.current_time.jd,
                            self.ephemeris['Time'].jd,
                            self.ephemeris['delta']) * units.au


Ephmeris = Body


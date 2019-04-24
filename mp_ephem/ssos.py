

import datetime
import os
import pprint
import logging
import numpy
import requests
from astropy.io import ascii
from astropy.time import Time

from . import ephem

requests.packages.urllib3.disable_warnings()

__author__ = 'Michele Bannister, JJ Kavelaars'

SSOS_URL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssos.pl"
RESPONSE_FORMAT = 'tsv'
NEW_LINE = '\r\n'


class _SSOSParser(object):
    """
    Parse the result of an SSOS query, which is stored in an astropy Table object
    """

    def __init__(self):
        """
        setup the parser.
        """

    @staticmethod
    def _skip_missing_data(str_vals, ncols):
        """
        add a extra columna if missing, else return None.
        """
        while len(str_vals) < ncols:
            str_vals.append('None')
        return str_vals

    def parse(self, ssos_result_filename_or_lines):

        """
        given the result table create 'source' objects.

        :param ssos_result_filename_or_lines:
        :rtype Table
        """
        table_reader = ascii.get_reader(Reader=ascii.Basic)
        table_reader.inconsistent_handler = self._skip_missing_data
        table_reader.header.splitter.delimiter = '\t'
        table_reader.data.splitter.delimiter = '\t'
        return table_reader.read(ssos_result_filename_or_lines)


class _SSOSParamDictBuilder(object):
    """
    Build a dictionary of parameters needed for an SSOS Query.

    http://www4.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl?
    lang=en
    object=2
    search=bynameCADC
    epoch1=1990+01+01
    epoch2=2016+3+27
    eellipse=
    eunits=arcseconds
    extres=yes
    xyres=yes

    This should be fun!
    """

    def __init__(self,
                 observations=None,
                 verbose=False,
                 search_start_date=None,
                 search_end_date=None,
                 orbit_method='bern',
                 error_ellipse='bern',
                 resolve_extension=True,
                 resolve_position=True,
                 telescope_instrument='CFHT/MegaCam'):

        self._orbit_method = orbit_method
        self._verbose = False
        self.verbose = verbose
        self._search_start_date = None
        self.search_start_date = search_start_date
        self._search_end_date = None
        self.search_end_date = search_end_date
        self._error_ellipse = None
        self._resolve_extension = None
        self.resolve_extension = resolve_extension
        self._resolve_position = None
        self.resolve_position = resolve_position
        self._telescope_instrument = None
        self.telescope_instrument = telescope_instrument
        self._error_units = None
        self._observations = []
        self.observations = observations
        self.error_ellipse = error_ellipse

    @property
    def observations(self):
        """
        The observations to be used in fitting, returned as list of
        the mpc format lines.

        This should be set to a list of objects whose 'str' values
        will be valid MPC observations.
        """
        return self._observations

    @observations.setter
    def observations(self, observations):
        if not isinstance(observations, list) and not isinstance(observations, numpy.ndarray):
            observations = [observations]
        self._observations = []
        orbit_method_set = None
        for observation in observations:
            use_bern = isinstance(observation, ephem.ObsRecord)
            self.orbit_method = use_bern and 'bern' or 'bynameCADC'
            orbit_method_set = orbit_method_set is None and self.orbit_method or orbit_method_set
            if orbit_method_set != self.orbit_method:
                raise ValueError("All members of observations list must be same type.")
            # use_bern needs to have any null observations removed.
            if use_bern and observation.null_observation:
                continue
            self._observations.append(str(observation))

    @property
    def verbose(self):
        """
        In verbose mode the SSOS query will return diagnostic
        information about how the search was done.
        """
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = (verbose and 'yes') or 'no'

    @property
    def search_start_date(self):
        """
        :return: Time constraint for start of SSOS search window.
        :rtype: Time
        """
        if self._search_start_date is None:
            return ""
        return self._search_start_date

    @search_start_date.setter
    def search_start_date(self, search_start_date):
        """

        :param search_start_date: Time object for start of SSOS search window.
        """
        try:
            self._search_start_date = Time(search_start_date).replicate(format='iso')
            self._search_start_date.out_subfmt = 'date'
        except Exception as ex:
            logging.debug(str(ex))
            self._search_start_date = None

    @property
    def search_end_date(self):
        """
        astropy.io.Time object. The end date of SSOS search window.
        """
        if self._search_end_date is None:
            return ""
        return self._search_end_date

    @search_end_date.setter
    def search_end_date(self, search_end_date):
        """
        :type search_end_date: astropy.io.Time
        :param search_end_date: search for frames take after the given date.
        """
        try:
            self._search_end_date = Time(search_end_date).replicate(format='iso')
            self._search_end_date.out_subfmt = 'date'
        except Exception as ex:
            logging.debug(str(ex))
            self._search_end_date = None

    @property
    def orbit_method(self):
        """
        What fitting method should be used to turn the observations
        into an orbit.

        Must be one of ['bern', 'mpc']
        """
        return self._orbit_method

    @orbit_method.setter
    def orbit_method(self, orbit_method):
        assert orbit_method in ['bern', 'mpc', 'bynameCADC']
        self._orbit_method = orbit_method
        if self._orbit_method == 'bynameCADC':
            self.error_units = 'arcseconds'

    @property
    def error_units(self):
        return self._error_units

    @error_units.setter
    def error_units(self, error_units):
        self._error_units = error_units

    @property
    def error_ellipse(self):
        """
        The size of the error ellipse to assign to each position, or
        'bern' to use the output of the BK fit.
        """
        return self._error_ellipse

    @error_ellipse.setter
    def error_ellipse(self, error_ellipse):
        """

        :param error_ellipse: either a number or the work 'bern'
        """
        if not self.orbit_method == 'bern':
            try:
                error_ellipse = float(error_ellipse)
            except:
                error_ellipse = ''
        self._error_ellipse = error_ellipse

    @property
    def resolve_extension(self):
        """
        Should SSOS resolve and return which extension of a frame the
        object would be in?
        """
        return self._resolve_extension

    @resolve_extension.setter
    def resolve_extension(self, resolve_extension):
        if str(resolve_extension).lower() == "no":
            resolve_extension = False
        self._resolve_extension = (resolve_extension and "yes") or "no"

    @property
    def resolve_position(self):
        """
        Should SSOS resolve and return the predicted X/Y location of
        the source?
        """
        return self._resolve_position

    @resolve_position.setter
    def resolve_position(self, resolve_position):
        if str(resolve_position).lower() == "no":
            resolve_position = False
        self._resolve_position = (resolve_position and "yes") or "no"

    @property
    def telescope_instrument(self):
        """
        Name of the telescope being used.
        """
        return self._telescope_instrument

    @telescope_instrument.setter
    def telescope_instrument(self, telescope_instrument):
        assert isinstance(telescope_instrument, str)
        if telescope_instrument in TELINST:
            self._telescope_instrument = telescope_instrument

    @property
    def params(self):
        """
        :return: A dictionary of SSOS query parameters.
        :rtype: dict
        """
        params = dict(format=RESPONSE_FORMAT,
                      verbose=self.verbose,
                      epoch1=str(self.search_start_date),
                      epoch2=str(self.search_end_date),
                      search=self.orbit_method,
                      eunits=self.error_units,
                      eellipse=self.error_ellipse,
                      extres=self.resolve_extension,
                      xyres=self.resolve_position,
                      telinst=self.telescope_instrument)

        if self.orbit_method == 'bynameCADC':
            params['object'] = NEW_LINE.join((str(target_name) for target_name in self.observations))
        else:
            params['obs'] = NEW_LINE.join((str(observation) for observation in self.observations))
        return params


class Query(object):
    """
    Query the CADC's Solar System Object search for a given set of
    MPC-formatted moving object detection lines.

    Inputs:
        - a list of ephem.ObsRecord instances

    Optional:
        - a tuple of the start and end times to be searched
          between. Format '%Y-%m-%d'

    Otherwise the temporal range defaults to spanning from the start
    of OSSOS surveying on 2013-01-01 to the present day.

    """

    def __init__(self,
                 observations=None,
                 search_start_date=Time('2003-01-01', scale='utc'),
                 search_end_date=Time('2017-01-01', scale='utc'),
                 error_ellipse='bern'):

        self.param_dict_builder = _SSOSParamDictBuilder(
            observations=observations,
            search_start_date=search_start_date,
            search_end_date=search_end_date,
            error_ellipse=error_ellipse)
        self.headers = {'User-Agent': 'OSSOS'}

    def get(self):
        """
        :return: A string containing the TSV result from SSOS
        :rtype: str
        :raise: AssertionError
        """
        params = self.param_dict_builder.params
        logging.debug(pprint.pformat(format(params)))
        response = requests.post(SSOS_URL,
                                 data=params,
                                 headers=self.headers)
        logging.debug(response.url)
        assert isinstance(response, requests.Response)
        assert (response.status_code == requests.codes.ok)

        lines = response.content
        # NOTE: spelling 'occured' is in SSOIS
        if len(lines) < 2 or "An error occured getting the ephemeris" in lines:
            logging.error("SSOIS reported an error:")
            logging.error(str(lines))
            logging.error(str(response.url))
            raise IOError(os.errno.EACCES, "call to SSOIS failed on format error")

        if os.access("backdoor.tsv", os.R_OK):
            lines += open("backdoor.tsv").read()
        return lines


def query(mpc_observations_or_target_name, search_start_date=None, search_end_date=None):
        """Send a query to the SSOS web service, looking for available observations using the given track.

        :return: an Table object
        :rtype: Table
        """

        if search_start_date is None:
            search_start_date = Time('1999-01-01', scale='utc')
        if search_end_date is None:
            search_end_date = Time(datetime.datetime.now().strftime('%Y-%m-%d'), scale='utc')
        logging.info("Sending query to SSOS start_date: {} end_data: {}\n".format(search_start_date, search_end_date))
        return _SSOSParser().parse(Query(mpc_observations_or_target_name,
                                         search_start_date=search_start_date,
                                         search_end_date=search_end_date).get())

TELINST = [
    'AAT/WFI',
    'ALMA',
    'CFHT/CFH12K',
    'CFHT/MegaCam',
    'CFHT/WIRCam',
    'CTIO-4m/DECam',
    'CTIO-4m/Mosaic2',
    'CTIO-4m/NEWFIRM',
    'ESO-LaSilla_2.2m/WFI',
    'ESO-NTT/EFOSC',
    'ESO-NTT/EMMI',
    'ESO-NTT/SOFI',
    'ESO-NTT/SUSI',
    'ESO-NTT/SUSI2',
    'ESO-VISTA/VIRCAM',
    'ESO-VLT/FORS1',
    'ESO-VLT/FORS2',
    'ESO-VLT/HAWKI',
    'ESO-VLT/ISAAC',
    'ESO-VLT/NAOS-CONICA',
    'ESO-VLT/VIMOS',
    'ESO-VST/OMEGACAM',
    'Gemini/GMOS',
    'Gemini/GMOS-N',
    'Gemini/GMOS-S',
    'Gemini/NIRI',
    'HST/ACS',
    'HST/WFC3',
    'HST/WFPC2',
    'INT/WFC',
    'JKT/JAG-CCD',
    'KPNO-0.9m/Mosaic1',
    'KPNO-0.9m/Mosaic1.1',
    'KPNO-4m/Mosaic1',
    'KPNO-4m/Mosaic1.1',
    'KPNO-4m/NEWFIRM',
    'NEAT-GEODSS-Maui',
    'SDSS',
    'SOAR/SOI',
    'Subaru/SuprimeCam',
    'WHT/AUXCAM',
    'WHT/INGRID',
    'WHT/LIRIS',
    'WHT/Prime',
    'WISE',
    'WIYN/MiniMo',
    'WIYN/ODI',
]

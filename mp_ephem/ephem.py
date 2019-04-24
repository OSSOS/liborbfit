# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import re
import struct
import logging
from astropy.coordinates import SkyCoord
from astropy import units
import numpy
from astropy.time import Time

__author__ = 'jjk, mtb55'

DEFAULT_OBSERVERS = ['M.P. Centre',
                     ]
DEFAULT_MEASURERS = ['I.C. You',
                     ]
DEFAULT_TELESCOPE = "OBSER+Size CCD"
DEFAULT_ASTROMETRIC_NETWORK = "UCAC4"

_KNOWN_OBSERVAOTRY_CODES = []
__PATH__ = os.path.dirname(__file__)
OBSERVATORIES_FILENAME = os.getenv('ORBIT_OBSERVATORIES', os.path.join(__PATH__, 'data', 'observatories.dat'))

# Build a list of observatories our observatories.dat file knows about, we will set the remaining ones to 500
observatories = open(OBSERVATORIES_FILENAME, mode='r', encoding='utf-8')
for line in observatories.readlines():
    if line.startswith('#'):
        continue
    _KNOWN_OBSERVAOTRY_CODES.append(line.split()[0])
observatories.close()

MPCNOTES = {"Note1": {" ": " ",
                      "": " ",
                      "*": "*",
                      "A": "earlier approximate position inferior",
                      "a": "sense of motion ambiguous",
                      "B": "bright sky/black or dark plate",
                      "b": "bad seeing",
                      "c": "crowded star field",
                      "D": "declination uncertain",
                      "d": "diffuse image",
                      "E": "at or near edge of plate",
                      "F": "faint image",
                      "f": "involved with emulsion or plate flaw",
                      "G": "poor guiding",
                      "g": "no guiding",
                      "H": "hand measurement of CCD image",
                      "I": "involved with star",
                      "i": "inkdot measured",
                      "J": "J2000.0 rereduction of previously-reported position",
                      "K": "stacked image",
                      "k": "stare-mode observation by scanning system",
                      "M": "measurement difficult",
                      "m": "image tracked on object motion",
                      "N": "near edge of plate, measurement uncertain",
                      "O": "image out of focus",
                      "o": "plate measured in one direction only",
                      "P": "position uncertain",
                      "p": "poor image",
                      "R": "right ascension uncertain",
                      "r": "poor distribution of reference stars",
                      "S": "poor sky",
                      "s": "streaked image",
                      "T": "time uncertain",
                      "t": "trailed image",
                      "U": "uncertain image",
                      "u": "unconfirmed image",
                      "V": "very faint image",
                      "W": "weak image",
                      "w": "weak solution"},
            "Note2": {" ": " ",
                      "": " ",
                      "P": "Photographic",
                      "e": "Encoder",
                      "C": "CCD",
                      "T": "Meridian or transit circle",
                      "M": "Micrometer",
                      "V": "'Roving Observer' observation",
                      "R": "Radar observation",
                      "S": "Satellite observation",
                      "c": "Corrected-without-republication CCD observation",
                      "E": "Occultation-derived observations",
                      "O": "Offset observations (used only for observations of natural satellites)",
                      "H": "Hipparcos geocentric observations",
                      "N": "Normal place",
                      "n": "Mini-normal place derived from averaging observations from video frames"},
            'PhotometryNote': {" ": " ",
                               "": " ",
                               "L": "Photometry uncertainty lacking",
                               "Y": "Photometry measured successfully",
                               "Z": "Photometry measurement failed."}}


class MinorPlanetNumber(object):
    PACKED_DESIGNATION = " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    def __init__(self, minor_planet_number):
        self._minor_planet_number = None
        if minor_planet_number is not None and len(minor_planet_number.strip(' ')) != 0:
            if len(minor_planet_number) == 5:
                if not minor_planet_number[0].isdigit():
                    base = 100000 * (MinorPlanetNumber.PACKED_DESIGNATION.index(minor_planet_number[0]))
                else:
                    base = 10000 * int(minor_planet_number[0])
                self._minor_planet_number = int(minor_planet_number[1:]) + base
            else:
                raise MPCFieldFormatError("minor_planet_number",
                                          "1 char or digit and 4 digits",
                                          minor_planet_number)

    def __format__(self, format_spec):
        return self.__str__()


    def __str__(self):
        if self._minor_planet_number is None:
            return ""
        if self._minor_planet_number > 99999:
            idx = int(self._minor_planet_number / 100000)
            minor_planet_number = MinorPlanetNumber.PACKED_DESIGNATION[idx]
            minor_planet_number += "{:04d}".format(self._minor_planet_number - idx * 100000)
        else:
            minor_planet_number = "{:05d}".format(self._minor_planet_number)
        return minor_planet_number

    def __int__(self):
        return self._minor_planet_number

    def __lt__(self, other):
        return int(self) < int(other)

    def __gt__(self, other):
        return int(self) > int(other)

    def __cmp__(self, other):
        return cmp(int(self), int(other))

    def __bool__(self):
        return self._minor_planet_number is not None

    def __nonzero__(self):
        return self.__bool__()


class NullObservation(object):
    NULL_OBSERVATION_CHARACTERS = ["!", "-", "#"]

    def __init__(self, null_observation=None, null_observation_character=None):
        """
        A boolean object that keeps track of True/False status via a set of magic characters.
        """
        if null_observation is not None and \
                isinstance(null_observation, str) and len(str(null_observation).strip(' ')) > 0 and \
                null_observation not in NullObservation.NULL_OBSERVATION_CHARACTERS:
            raise MPCFieldFormatError("null_observation",
                                      "one of " + str(NullObservation.NULL_OBSERVATION_CHARACTERS),
                                      null_observation)
        if null_observation_character is None:
            null_observation_character = NullObservation.NULL_OBSERVATION_CHARACTERS[0]
        self.null_observation_character = null_observation_character

        self._null_observation = None
        if isinstance(null_observation, str):
            self._null_observation = null_observation in NullObservation.NULL_OBSERVATION_CHARACTERS
        elif isinstance(null_observation, bool):
            self._null_observation = null_observation
        else:
            self._null_observation = False

    def __format__(self, format_spec):
        return self.__str__()

    def __str__(self):
        return self._null_observation and self.null_observation_character or " "

    def __bool__(self):
        return self._null_observation

    def __nonzero__(self):
        return self.__bool__()


class MPCFormatError(Exception):
    """Base class for errors in MPC formatting."""


class TNOdbFlags(object):
    """
    The OSSOS/CFEPS database has a 'flag' field that indicates OSSOS specific issues associated with an
    MPC formatted line in the database.
    """

    def __init__(self, flags):
        if not re.match("\d{12}", flags):
            raise ValueError("illegal flag string: {}".format(flags))
        self.__flags = flags

    def __format__(self, format_spec):
        return self.__str__()

    def __str__(self):
        return self.__flags

    @property
    def is_discovery(self):
        """
        Is this observation part of the discovery triplet?  bit 1
        :return: bool
        """
        return self.__flags[0] == 1

    @is_discovery.setter
    def is_discovery(self, is_discovery):
        self.__flags[0] == bool(is_discovery) and "1" or "0"

    @property
    def is_secret(self):
        """
        Is this observation secret? bit 2
        :return: bool
        """
        return self.__flags[1] == 1


class MPCFieldFormatError(MPCFormatError):
    def __init__(self, field, requirement, actual):
        super(MPCFieldFormatError, self).__init__(
            "Field %s: %s; but was %s" % (field, requirement, actual))


def format_ra_dec(ra_deg, dec_deg):
    """
    Converts RA and DEC values from degrees into the formatting required
    by the Minor Planet Center:

    Formats:
      RA: 'HH MM SS.ddd'
      DEC: 'sDD MM SS.dd' (with 's' being the sign)

    (From: http://www.minorplanetcenter.net/iau/info/OpticalObs.html)

    Args:
      ra_deg: float
        Right ascension in degrees
      dec_deg: float
        Declination in degrees

    Returns:
      formatted_ra: str
      formatted_dec: str
    """
    coords = SkyCoord(ra=ra_deg, dec=dec_deg,
                      unit=(units.degree, units.degree))

    # decimal=False results in using sexagesimal form
    formatted_ra = coords.ra.to_string(unit=units.hour, decimal=False,
                                       sep=" ", precision=3, alwayssign=False,
                                       pad=True)

    formatted_dec = coords.dec.to_string(unit=units.degree, decimal=False,
                                         sep=" ", precision=2, alwayssign=True,
                                         pad=True)

    return formatted_ra, formatted_dec


class MPCNote(object):
    """
    Alphabetic note shown with some of the observations. Non-alphabetic codes are used to differentiate between
    different programs at the same site and such codes will be defined in the headings for the individual
    observatories in the Minor Planet Circulars.
    """

    def __init__(self, code="C", note_type="Note2"):
        self._note_type = None
        self._code = None
        self.note_type = note_type
        self.code = code

    @property
    def note_type(self):
        """
        Note 1 or 2 from an MPC line.
        """
        return self._note_type

    @note_type.setter
    def note_type(self, note_type):
        if note_type not in list(MPCNOTES.keys()):
            raise ValueError("Invalid note_type: expected one of %s got %s" % (str(list(MPCNOTES.keys())), note_type))
        self._note_type = note_type

    @property
    def code(self):
        """
        The MPC note the denotes the type of detector system used
        """
        return self._code

    @code.setter
    def code(self, code):
        """

        :type code: str
        :param code: an MPC Note code. Either from the allow dictionary or 0-9
        """
        if code is None:
            _code = " "
        else:
            _code = str(code).strip()

        if _code.isdigit():
            if self.note_type != 'Note1':
                raise MPCFieldFormatError(self.note_type,
                                          "Must be a character",
                                          _code)
            if not 0 < int(_code) < 10:
                raise MPCFieldFormatError(self.note_type,
                                          "numeric value must be between 0 and 9",
                                          _code)
        else:
            if len(_code) > 1:
                raise MPCFieldFormatError(self.note_type,
                                          "must be 0 or 1 characters",
                                          _code)
            if _code not in MPCNOTES[self.note_type]:
                logging.warning("Unknown note value: {}".format(_code))
        self._code = _code

    def __format__(self, format_spec):
        return self.__str__()

    def __str__(self):
        return str(self.code)

    @property
    def long(self):
        return MPCNOTES[self.note_type][self.code]


class Discovery(object):
    """
    Holds the discovery flag for an MPC ObsRecord Line
    """

    def __init__(self, is_discovery=False):
        self._is_discovery = False
        self._is_initial_discovery = False
        self.is_discovery = is_discovery
        self.is_initial_discovery = is_discovery

    def set_from_mpc_line(self, mpc_line):
        """
        Given an MPC line set the discovery object
        """
        mpc_line = str(mpc_line)
        if len(mpc_line) < 56:
            raise MPCFieldFormatError("mpc_line",
                                      "is too short",
                                      mpc_line)
        self.is_discovery = mpc_line[12]

    @property
    def is_initial_discovery(self):
        return self._is_initial_discovery

    @property
    def is_discovery(self):
        return self._is_discovery

    @is_discovery.setter
    def is_discovery(self, is_discovery):
        if is_discovery not in ['*', '&', ' ', '', True, False, None]:
            raise MPCFieldFormatError("discovery",
                                      "must be one of '',' ','&', '*',True, False. ",
                                      is_discovery)
        self._is_discovery = is_discovery in ['*', '&', True] and True or False

    @is_initial_discovery.setter
    def is_initial_discovery(self, is_discovery):
        """
        Is this MPC line the initial discovery line?
        @param is_discovery: the code for the discovery setting "*" or True or False
        """
        self._is_initial_discovery = (is_discovery in ["*", True] and True) or False

    def __str__(self):
        if self.is_initial_discovery:
            return "*"
        if self.is_discovery:
            return "&"
        return " "

    def __bool__(self):
        return self.is_discovery

    def __nonzero__(self):
        return self.__bool__()


def compute_precision(coord):
    """
    Returns the number of digits after the last '.' in a given number or string.

    """
    coord = str(coord).strip(' ')
    idx = coord.rfind('.')
    precision = 0
    if idx > 0:
        precision = len(coord) - idx - 1
    return precision

def get_date(date_string):
    """
    Given an MPC formatted time string return a Time object.

    :rtype : Time
    :param date_string: a string in MPC date format
    """
    _date_precision = compute_precision(date_string)
    return Time(date_string, format='mpc', scale='utc', precision=_date_precision)


class ObserverLocation(object):
    """
    Store the geocentric location of the observer
    """

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class ObsRecord(object):
    """
    An observation of an object, nominally generated by reading an MPC formatted file.
    """

    def __init__(self,
                 null_observation=False,
                 minor_planet_number=None,
                 provisional_name=None,
                 discovery=False,
                 note1=None,
                 note2=None,
                 date="2000 01 01.000001",
                 ra="00 00 00.000",
                 dec="+00 00 00.00",
                 mag=-1,
                 band='r',
                 observatory_code=568,
                 comment=None,
                 mag_err=-1,
                 xpos=None,
                 ypos=None,
                 frame=None,
                 plate_uncertainty=None,
                 astrometric_level=0):
        """

        :param provisional_name:
        :param discovery:
        :param note1:
        :param note2:
        :param date:
        :param ra:
        :param dec:
        :param mag:
        :param band:
        :param observatory_code:
        :param comment: A comment about this observation, not sent
        :param mag_err:
        :param xpos:
        :param ypos:
        :param frame:
        :param plate_uncertainty:
        :param null_observation:

        """
        self._null_observation = False
        self._minor_planet_number = None
        self._provisional_name = None

        self.null_observation = null_observation
        self.minor_planet_number = minor_planet_number
        self.provisional_name = provisional_name

        self._discovery = None
        self.discovery = discovery
        self._note1 = None
        self.note1 = note1
        self._note2 = None
        self.note2 = note2
        self._date = None
        self._date_precision = None
        self.date = date
        self._coordinate = None
        self._mag = None
        self._mag_err = None
        self._mag_precision = 1
        self._ra_precision = 3
        self._dec_precision = 2
        self.coordinate = (ra, dec)
        self.mag = mag
        self.mag_err = mag_err
        self._band = None
        self.band = band
        self._observatory_code = None
        self.observatory_code = observatory_code
        self._comment = None
        self.location = None
        self.comment = OSSOSComment(version="O", frame=frame,
                                    source_name=provisional_name,
                                    photometry_note="",
                                    mpc_note=str(self.note1),
                                    x=xpos,
                                    y=ypos,
                                    plate_uncertainty=plate_uncertainty,
                                    astrometric_level=astrometric_level,
                                    magnitude=mag,
                                    mag_uncertainty=mag_err,
                                    comment=comment)

    def __eq__(self, other):
        return self.to_string() == other.to_string()

    def __ne__(self, other):
        return self.to_string() != other.to_string()

    def __le__(self, other):
        return self.date <= other.date

    def __lt__(self, other):
        return self.date < other.date

    def __ge__(self, other):
        return self.date >= other.date

    def __gt__(self, other):
        return self.date > other.date

    @classmethod
    def from_ted(cls, ted):
        """
        Turn a ted line into an MPC formatted line. (Marc Buie's format)

        Example line: 2011 04  28.29389  18 32 37.260  -21 13 49.86  26.3R NI100     304
        :param ted: a TED formatted minor planet observations
        :return: ObsRecord
        """
        ted = ted.strip()
        if len(ted) != len('2011 04  28.29389  18 32 37.260  -21 13 49.86  26.3R NI100     304'):
            raise ValueError("Incorrectly formatted line for ted to mpc conversion.")
        line_order = ["date", "ra", "dec", "mag", "filter", "provisional_name", "observatory_code"]
        parts = {"date": "2011 04  28.29389  ",
                 "ra": "18 32 37.260  ",
                 "dec": "-21 13 49.86  ",
                 "mag": "26.3",
                 "filter": "R ",
                 "provisional_name": "NI100     ",
                 "observatory_code": "304"}
        end_pos = 0
        args = {}
        for part in line_order:
            start_pos = end_pos
            end_pos = start_pos + len(parts[part])
            args[part] = ted[start_pos:end_pos]
        return ObsRecord(provisional_name=args['provisional_name'],
                         date=args['date'].replace("  ", " "),
                         ra=args['ra'],
                         dec=args['dec'],
                         mag=args['mag'],
                         band=args['filter'],
                         observatory_code=args['observatory_code'])

    @classmethod
    def from_string(cls, input_line):
        """
        Given an MPC formatted line, returns an MPC ObsRecord object.

        This method attempts four different format.

        If line is less than 80 characters then assumes in 'ted' format,

        If line equal/greater than 80 tries some fixed structure formats

        mpc_format = '0s5s7s1s1s1s17s12s12s9x5s1s6x3s'
        ossos_format1 = '1s4s7s1s1s1s17s12s12s9x5s1s6x3s'
        ossos_format2' = '1s0s11s1s1s1s17s12s12s9x5s1s6x3s'

        if those fail then tries Alex Parker's .ast format.

        :param input_line: a line in the one-line roving observer format
        :type input_line: str
        """
        struct_formats = {'mpc_format': '0s5s7s1s1s1s17s12s12s9x5s1s6x3s',
                          'ossos_format1': '1s4s7s1s1s1s17s12s12s9x5s1s6x3s',
                          'ossos_format2': '1s0s11s1s1s1s17s12s12s9x5s1s6x3s'}

        mpc_line = input_line.strip('\n')
        logging.debug("Trying to create MPC record from:\n{}".format(mpc_line))
        if len(mpc_line) > 0 and mpc_line[0] == '#':
            return MPCComment.from_string(mpc_line[1:])

        if mpc_line[31:34] == ' 1 ':
            # This is a spacecraft telemetry line for the previous Observation.
            x = float(mpc_line[34] + mpc_line[35:45].strip())
            y = float(mpc_line[46] + mpc_line[47:57].strip())
            z = float(mpc_line[58] + mpc_line[59:69].strip())
            return ObserverLocation(x, y, z)

        comment = mpc_line[81:]
        mpc_line = mpc_line[0:80]
        if len(mpc_line) != 80 and len(mpc_line) > 0:
            logging.info("{}".format(mpc_line))
            logging.info("mpc line is only {} chars long, trying .ted format".format(len(mpc_line)))
            try:
                return cls.from_ted(mpc_line)
            except Exception as ex:
                logging.debug(type(ex))
                logging.debug(str(ex))

        obsrec = None
        for format_name in struct_formats:
            try:
                logging.debug("Trying to parse with {}".format(format_name))
                args = [str(arg, 'utf-8').strip() for arg in struct.unpack(struct_formats[format_name],
                                                                           bytes(mpc_line, 'utf-8'))]
                obsrec = cls(*args)
                break
            except Exception as ex:
                logging.debug("Failed: {}".format(ex))
                obsrec = None

        if obsrec is None:
            logging.debug("Trying AP's .ast format")
            # try converting using Alex Parker's .ast format:
            # 2456477.78468 18:39:07.298 -20:40:17.53 0.2 304
            try:
                _parts = input_line.split(' ')
                args = {"date": Time(float(_parts[0]), scale='utc', format='jd').mpc,
                        "discovery": False,
                        "ra": _parts[1].replace(":", " "),
                        "dec": _parts[2].replace(":", " "),
                        "plate_uncertainty": _parts[3],
                        "observatory_code": _parts[4]}
                obsrec = cls(**args)
                return obsrec
            except Exception as ex:
                logging.debug("Failed to parse as AP format: {}".format(str(ex)))

        if obsrec is None:
            logging.debug("Trying Simon Porter format")
            try:
                _parts = input_line.split()
                year = int(_parts[1])
                month = int(_parts[2])
                day = float(_parts[3])
                obs_date = Time("{:4d} {:02d} {:08.5f}".format(year, month, day), format='mpc', precision=5)
                args = {"date": obs_date,
                        "provisional_name": _parts[0],
                        'discovery': False,
                        'ra': "{} {} {}".format(_parts[7], _parts[8], _parts[9]),
                        'dec': "{} {} {}".format(_parts[10], _parts[11], _parts[12]),
                        'mag': float(_parts[13]),
                        'mag_err': float(_parts[14]),
                        'observatory_code': 500,
                        'comment': " ".join(_parts[15:])}
                obsrec = cls(**args)
                obsrec.location = ObserverLocation(_parts[4], _parts[5], _parts[6])
                return obsrec
            except Exception as ex:
                logging.debug("Failed to parse as Simon Porter format: {}".format(str(ex)))

        if obsrec is None or not obsrec:
            if mpc_line is not None and len(mpc_line) > 0:
                logging.warning("Failed to parse line: {}".format(mpc_line))
            return obsrec

        obsrec.comment = MPCComment.from_string(comment)
        if isinstance(obsrec.comment, OSSOSComment) and obsrec.comment.source_name is None:
            obsrec.comment.source_name = obsrec.provisional_name
        # Check if there are TNOdb style flag lines.
        if hasattr(obsrec.comment, 'flags') and obsrec.comment.flags is not None:
            try:
                if obsrec.comment.flags[0] == '1':
                    obsrec.discovery.is_discovery = True
            except Exception as ex:
                logging.warning(str(ex))
                logging.warning("Invalid flag string in comment: {}".format(obsrec.comment.flags))

        return obsrec

    def to_string(self):
        as_string = str(self)
        if self.comment is not None and str(self.comment) != "":
            as_string += " " + str(self.comment)
        return as_string

    def __str__(self):
        """
        Writes out data about accepted objects in the Minor Planet Center's 'Minor Planets'
        format as specified here:
        http://www.minorplanetcenter.net/iau/info/OpticalObs.html
        """
        # MOP/OSSOS allows the provisional name to take up the full space allocated to the MinorPlanetNumber AND
        # the provisional name.

        if not self.minor_planet_number:
            if len(self.provisional_name) > 7:
                mpc_str = "{:1.1s}{:<11.11s}".format(self.null_observation,
                                                     self.provisional_name)
            else:
                mpc_str = "{:1.1s}{:4.4s}{:<7.7s}".format(self.null_observation,
                                                          " " * 4,
                                                          self.provisional_name)
        else:
            mpc_str = "{:5.5s}{:<7.7s}".format(self.minor_planet_number, self.provisional_name)

        mpc_str += str(self.discovery)
        mpc_str += '{0:1s}{1:1s}'.format(str(self.note1), str(self.note2))
        mpc_str += '{0:<17s}'.format(str(self.date))
        mpc_str += '{0:<12s}{1:<12s}'.format(str(self.ra), str(self.dec))
        mpc_str += 9 * " "
        mag_format = '{0:<5.' + str(self._mag_precision) + 'f}{1:1s}'
        band = self.band is None and " " or self.band
        mag_str = (self.mag is None and 6 * " ") or mag_format.format(self.mag, band)
        if len(mag_str) != 6:
            raise MPCFieldFormatError("mag",
                                      "length of mag string should be exactly 6 characters, got->",
                                      mag_str)
        mpc_str += mag_str
        mpc_str += 6 * " "
        mpc_str += "%3s" % self.observatory_code
        return mpc_str

    def to_tnodb(self):
        """
        provide string representation of observation in a format used for OSSOS database input.
        """

        # O indicates OSSOS survey
        if not isinstance(self.comment, OSSOSComment):
            logging.warn("Non OSSOS comment:{}".format(self.comment))

        comment_line = "#" + str(self.comment).rstrip('\n')

        if self.mag == -1:  # write no mag and no filter for where photometry couldn't be measured
            self.mag = None
        else:
            # set mag precision back to 0.1 mags regardless of how good it actually is
            self._mag_precision = 1

        # set the null observation character to the tnodb value
        self.null_observation.null_observation_character = "-"
        mpc_observation = str(self)

        return comment_line + '\n' + mpc_observation

    def to_mpc(self):
        self.null_observation.null_observation_character = "#"
        return str(self)

    @property
    def null_observation(self):
        return self._null_observation

    @null_observation.setter
    def null_observation(self, null_observation=False):
        """
        :param null_observation: is this a null observation marker True/False
        """
        self._null_observation = NullObservation(null_observation)

    @property
    def provisional_name(self):
        return self._provisional_name

    @provisional_name.setter
    def provisional_name(self, provisional_name=None):
        if provisional_name is None:
            provisional_name = " " * 7
        else:
            provisional_name = provisional_name.strip()
            # if not provisional_name[0].isalpha():
            # logging.warning("Provisional Name should not be a number: {}".format(provisional_name))
            # if not len(provisional_name) <= 7:
            #     logging.warning("Provisional Name too long {}".format(provisional_name))
        self._provisional_name = provisional_name

    @property
    def minor_planet_number(self):
        """
        :return: minor_planet_number for object associated with this observation.
        """
        return self._minor_planet_number

    @minor_planet_number.setter
    def minor_planet_number(self, minor_planet_number):
        self._minor_planet_number = MinorPlanetNumber(minor_planet_number)

    @property
    def discovery(self):
        """
        Is this a discovery observation?

        :return True/False
        """
        return self._discovery

    @discovery.setter
    def discovery(self, is_discovery):
        """

        :type is_discovery: bool
        :param is_discovery: indicates if observation was a discovery
        """
        self._discovery = Discovery(is_discovery=is_discovery)

    @property
    def note1(self):
        return self._note1

    @note1.setter
    def note1(self, note1):
        self._note1 = MPCNote(code=note1, note_type="Note1")

    @property
    def note2(self):
        return self._note2

    @note2.setter
    def note2(self, code):
        self._note2 = MPCNote(code=code, note_type="Note2")

    @property
    def date(self):
        return self._date

    @date.setter
    def date(self, date_str):
        self._date_precision = compute_precision(date_str)
        logging.debug("Setting precision to: {}".format(self._date_precision))
        try:
            self._date = Time(date_str, format='mpc', scale='utc', precision=self._date_precision)
        except Exception as ex:
            logging.warning(str(ex))
            raise MPCFieldFormatError("ObsRecord Date",
                                      "does not match expected format",
                                      date_str)

    @property
    def ra(self):
        return self.coordinate.ra.to_string(unit=units.hour, decimal=False,
                                            sep=" ", precision=self._ra_precision, alwayssign=False,
                                            pad=True)

    @property
    def dec(self):
        return self.coordinate.dec.to_string(unit=units.degree, decimal=False,
                                             sep=" ", precision=self._dec_precision, alwayssign=True,
                                             pad=True)

    @property
    def comment(self):
        """

        :return: the comment
        :rtype: OSSOSComment
        """
        return self._comment

    @comment.setter
    def comment(self, comment):
        if comment is None:
            self._comment = ""
        else:
            self._comment = comment

    @property
    def coordinate(self):
        return self._coordinate

    @coordinate.setter
    def coordinate(self, coord_pair):
        """

        :param coord_pair: RA/DEC pair [as a tuple or single string]
        """

        if type(coord_pair) in [list, tuple] and len(coord_pair) == 2:
            val1 = coord_pair[0]
            val2 = coord_pair[1]
        else:
            raise MPCFieldFormatError("RA/DEC",
                                      "Expected a pair of coordinates got: ",
                                      coord_pair)

        self._ra_precision = 3
        self._dec_precision = 2

        # First try using just the values as provided, maybe they have units
        try:
            self._coordinate = SkyCoord(val1, val2)
            return
        except Exception as ex:
            logging.debug("Failed to turn {} {} into SkyCoord".format(val1, val2))
            logging.debug(str(ex))
            pass

        # Now try treating them as floats in degrees.
        try:
            ra = float(val1)
            dec = float(val2)
            self._coordinate = SkyCoord(ra, dec, unit=(units.degree, units.degree))
            return
        except Exception as ex:
            logging.debug("Failed to create coordinate using RA/DEC as degrees: {}/{}".format(val1, val2))
            logging.debug(str(ex))
            pass

        # Maybe they are strings and then assume HH:MM:SS dd:mm:ss
        try:
            self._ra_precision = compute_precision(val1)
            self._ra_precision = self._ra_precision < 3 and self._ra_precision or 3
            self._dec_precision = compute_precision(val2)
            self._dec_precision = self._dec_precision < 2 and self._dec_precision or 2
            self._coordinate = SkyCoord(val1, val2, unit=(units.hour, units.degree))
            return
        except Exception as ex:
            logging.debug(type(ex))
            logging.debug(str(ex))
            raise MPCFieldFormatError("coord_pair",
                                      "must be [ra_deg, dec_deg] or HH MM SS.S[+-]dd mm ss.ss",
                                      coord_pair)

    @property
    def mag(self):
        return self._mag

    @mag.setter
    def mag(self, mag):
        if mag is None or len(str(str(mag).strip(' '))) == 0 or float(mag) < 0:
            self._mag_precision = 0
            self._mag = None
        else:
            self._mag = float(mag)
            self._mag_precision = min(1, compute_precision(str(mag)))

    @property
    def mag_err(self):
        return self._mag_err

    @mag_err.setter
    def mag_err(self, mag_err):
        if mag_err is None or len(str(mag_err).strip('')) == 0 or self.mag is None:
            self._mag_err = None
        else:
            self._mag_err = mag_err

    @property
    def band(self):
        return self._band

    @band.setter
    def band(self, band):
        band = str(band.strip(' '))
        self._band = (len(band) > 0 and str(band)[0]) or None

    @property
    def observatory_code(self):
        return self._observatory_code

    @observatory_code.setter
    def observatory_code(self, observatory_code):
        observatory_code = str(observatory_code)
        if observatory_code not in _KNOWN_OBSERVAOTRY_CODES:
            observatory_code = "500"
        if not len(observatory_code) <= 4:
            raise MPCFieldFormatError("Observatory code",
                                      "must be 3 characters or less",
                                      observatory_code)
        self._observatory_code = str(observatory_code)


# Create an Observation class that matches the ObsRecord class.  The new class name is ObsRecord (for disambiguation)
# but the old name is kept available for backward compatibility.
Observation = ObsRecord


class MPCComment(object):
    """
    A generic class for all comment strings.. try and figure out which one to use.
    """

    def __init__(self, comment):
        self.comment = str(comment)

    def __format__(self, format_spec):
        return self.__str__()

    def __str__(self):
        return str(self.comment)

    @classmethod
    def from_string(cls, line):
        comment = line
        logging.debug('Here is the comment. ' + comment)
        for func in [TNOdbComment.from_string,
                     OSSOSComment.from_string,
                     RealOSSOSComment.from_string,
                     CFEPSComment.from_string,
                     MPCComment]:
            try:
                comment = func(line)
            except ValueError as verr:
                logging.debug(str(func))
                logging.debug(str(verr))
                continue
            break
        return comment

    def to_string(self):
        as_string = str(self)
        if self.comment is not None and str(self.comment) != "":
            as_string += " " + str(self.comment)
        return as_string


class OSSOSComment(object):
    """
    Parses an OSSOS observation's metadata into a format that can be stored in the 
    an ObsRecord.comment and written out in the same MPC line.

    Specification: '1s1x12s1x11s1x2s1x7s1x7s1x4s1x1s1x5s1x4s1x'
    """

    def __init__(self, version, frame, source_name, photometry_note, mpc_note, x, y,
                 plate_uncertainty=0.2,
                 astrometric_level=0,
                 magnitude=None,
                 mag_uncertainty=None,
                 comment=None):
        self.version = version
        self._frame = None
        self.frame = frame
        self.source_name = source_name
        self._photometry_note = None
        self.photometry_note = photometry_note
        self.mpc_note = mpc_note
        self._x = None
        self.x = x
        self._y = None
        self.y = y
        self._mag = None
        self.mag = magnitude
        self._mag_uncertainty = None
        self.mag_uncertainty = mag_uncertainty
        self._plate_uncertainty = None
        self._astrometric_level = 0
        try:
            self.plate_uncertainty = plate_uncertainty
            self.astrometric_level = astrometric_level
        except:
            self.plate_uncertainty = 0.2
            self.mag = plate_uncertainty
            self.mag_uncertainty = astrometric_level
        self._comment = ""
        self.comment = comment
        self.flags = None

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return str(self) != str(other)

    def __le__(self, other):
        raise NotImplemented

    def __ge__(self, other):
        raise NotImplemented

    @classmethod
    def from_string(cls, comment):
        """
        Build an MPC Comment from a string.
        """
        logging.debug("Attempting to build an OSSOS comment from {}".format(comment))
        if comment is None or len(comment) == 0:
            return str("")
        if comment[0] == "#":
            comment = comment[1:]
        values = comment.split('%')
        comment_string = ""
        if len(values) > 1:
            comment_string = values[1].lstrip(' ')
        # O 1631355p21 O13AE2O     Z  1632.20 1102.70 0.21 3 ----- ---- % Apcor failure.
        ossos_comment_format = '1s1x12s1x11s1x1s1s1x7s1x7s1x4s1x1s1x5s1x4s1x'
        old_ossos_comment_format = '1s1x10s1x11s1x1s1s1x7s1x7s1x4s1x1s1x5s1x4s1x'
        bad_ossos_comment_format = '1s3x10s1x11s1x1s1s1x7s1x7s1x4s1x1s1x5s1x4s1x'
        orig_ossos_comment_format = '1s1x10s1x8s1x1s1s1x6s1x6s1x5s1x4s1x4s1x'

        for struct_ in [ossos_comment_format, old_ossos_comment_format, bad_ossos_comment_format,
                        orig_ossos_comment_format]:
            try:
                logging.debug("Parsing using structure: {}".format(struct_))
                args = [str(arg, 'utf-8') for arg in struct.unpack(struct_, bytes(values[0], 'utf-8'))]
                retval = cls(*args)
                retval.comment = values[1]
                return retval
            except Exception as e:
                logging.debug(str(e))
                logging.debug("OSSOS Fixed Format Failed.")

        logging.debug("Trying space separated version")
        values = values[0].split()
        # O 1645236p27 L3UV Y 106.35 4301.85 0.20 0 24.31 0.15
        try:
            if values[0] != 'O' or len(values) < 6:
                # this is NOT and OSSOS style comment string
                raise ValueError("Can't parse non-OSSOS style comment: {}".format(comment))
            # first build a comment based on the required fields.
            retval = cls(version="O",
                         frame=values[1].strip(),
                         source_name=values[2].strip(),
                         photometry_note=values[3][0].strip(),
                         mpc_note=values[3][1:].strip(),
                         x=values[4].strip(),
                         y=values[5].strip(),
                         plate_uncertainty=len(values) > 6 and float(values[6]) or None,
                         astrometric_level=len(values) > 7 and int(values[7]) or None,
                         magnitude=len(values) > 8 and float(values[8]) or None,
                         mag_uncertainty=len(values) > 9 and float(values[9]) or None,
                         comment=comment_string.strip())
        except Exception as e:
            logging.debug(values)
            logging.debug(comment)
            logging.debug(str(e))
            raise e

        retval.version = values[0]
        logging.debug("length of values: {}".format(len(values)))
        logging.debug("Values: {}".format(str(values)))
        # the format of the last section evolved during the survey, but the following flags should handle this.
        if len(values) == 7:
            retval.plate_uncertainty = values[6]
        elif len(values) == 8:
            retval.plate_uncertainty = values[6]
            retval.astrometric_level = values[7]
        elif len(values) == 9:  # This is the old format where mag was in-between X/Y and uncertainty in X/Y
            retval.mag = values[6]
            retval.mag_uncertainty = values[7]
            retval.plate_uncertainty = values[8]
        elif len(values) == 10:
            if float(values[8]) < 1:  # if there are 9 values then the new astrometric level value is set after mag
                retval.plate_uncertainty = values[8]
                retval.astrometric_level = values[9]
                logging.debug('here now')
                retval.mag = values[6]
                retval.mag_uncertainty = values[7]
        logging.debug("DONE.")
        return retval

    @property
    def frame(self):
        return self._frame

    @frame.setter
    def frame(self, frame):
        if frame is not None:
            self._frame = "{}".format(frame).strip()
        else:
            self._frame = frame

    @property
    def mag(self):
        return self._mag

    @mag.setter
    def mag(self, mag):
        try:
            self._mag = float(mag)
            self.photometry_note = "Y"
            if not 15 < self._mag < 30:
                raise ValueError("Magnitude out of reasonable range:  15 < mag < 30")
        except:
            self.photometry_note = "Z"
            self._mag = None

    @property
    def mag_uncertainty(self):
        return self._mag_uncertainty

    @mag_uncertainty.setter
    def mag_uncertainty(self, mag_uncertainty):
        try:
            self._mag_uncertainty = float(mag_uncertainty)
            if not 0 < self._mag_uncertainty < 1.0:
                raise ValueError("mag uncertainty must be in range 0 to 1")
        except Exception as e:
            logging.debug("Failed trying to convert merr ({}) to float. Using default.".format(mag_uncertainty))
            logging.debug(str(e))
            if str(self.mag).isdigit():
                self.photometry_note = "L"
            else:
                self.photometry_note = "Z"
            self._mag_uncertainty = None

    @property
    def photometry_note(self):
        return self._photometry_note

    @property
    def astrometric_level(self):
        return self._astrometric_level

    @astrometric_level.setter
    def astrometric_level(self, astrometric_level):
        try:
            astrometric_level = int(astrometric_level)
        except:
            astrometric_level = 0
        if not -1 < astrometric_level < 10:
            raise ValueError("Astrometric level must be integer between 0 and 9.")
        self._astrometric_level = astrometric_level

    @photometry_note.setter
    def photometry_note(self, photometry_note):
        self._photometry_note = str(photometry_note)

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        try:
            self._x = float(x)
        except Exception:
            self._x = None

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        try:
            self._y = float(y)
        except Exception:
            self._y = None

    @property
    def plate_uncertainty(self):
        return self._plate_uncertainty

    @plate_uncertainty.setter
    def plate_uncertainty(self, plate_uncertainty):
        try:
            self._plate_uncertainty = float(plate_uncertainty)
        except Exception:
            self._plate_uncertainty = 0.2
        if not 0 < self._plate_uncertainty < 1.:
            raise ValueError("Plate uncertainty must be between 0 and 1. (in arc-seconds)")

    @property
    def comment(self):
        return self._comment

    @comment.setter
    def comment(self, comment):
        if comment is not None:
            try:
                self._comment = str(comment.strip())
            except Exception:
                self._comment = ''
        else:
            self._comment = ''

    @staticmethod
    def to_str(frmt, value, default="", sep=" "):
        try:
            if value is None:
                raise ValueError("Don't print None.")
            return sep + frmt.format(value)
        except Exception:
            return sep + default

    def __format__(self, format_spec):
        return self.__str__()

    def __str__(self):
        """
        Format comment as required for storing OSSOS metadata
        """
        if self.version == "T":
            return self.comment

        if self.version == "L":
            return "{:1s} {:10s} {}".format(self.version, self.frame, self.comment)

        comm = '{:1s}'.format(self.version)
        comm += self.to_str("{:>12.12s}", self.frame, "-" * 12)
        comm += self.to_str("{:<11.11s}", self.source_name, "-" * 11)
        comm += self.to_str("{:2.2s}", self.photometry_note + self.mpc_note, "--")
        comm += self.to_str("{:>7.2f}", self.x, "-" * 7)
        comm += self.to_str("{:>7.2f}", self.y, "-" * 7)
        comm += self.to_str('{:4.2f}', self.plate_uncertainty, "-" * 4)
        comm += self.to_str('{:1d}', self.astrometric_level, "-")
        comm += self.to_str('{:5.2f}', self.mag, "-" * 5)
        comm += self.to_str('{:4.2f}', self.mag_uncertainty, "-" * 4)
        comm += ' % {}'.format(self.comment)

        return comm


class MPCWriter(object):
    """
    Writes out data about accepted objects in the Minor Planet Center's
    format as specified here:
    http://www.minorplanetcenter.net/iau/info/OpticalObs.html

    Note that we assume objects fall under the Minor Planet category.

    Format reproduced below for convenience:

        Columns     Format   Use
        1 -  5        A5     Minor planet number
        6 - 12        A7     Provisional or temporary designation
        13            A1     Discovery asterisk
        14            A1     Note 1
        15            A1     Note 2
        16 - 32      A17     Date of observation
        33 - 44      A12     Observed RA (J2000.0)
        45 - 56      A12     Observed Decl. (J2000.0)
        57 - 65       9X     Must be blank
        66 - 71    F5.2,A1   Observed magnitude and band
                               (or nuclear/total flag for comets)
        72 - 77       6X     Must be blank
        78 - 80       A3     Observatory code
    """

    def __init__(self, file_handle, auto_flush=True, include_comments=True,
                 auto_discovery=True, formatter=None):
        self.filehandle = file_handle
        self.auto_flush = auto_flush
        self.include_comments = include_comments

        # Holds observations that have not yet been flushed
        self.buffer = {}
        self._written_mpc_observations = []

        self.auto_discovery = auto_discovery
        self._discovery_written = False
        if formatter is None:
            if self.include_comments:
                self.formatter = ObsRecord.to_string
            else:
                self.formatter = ObsRecord.__str__
        else:
            self.formatter = formatter

    def get_filename(self):
        return self.filehandle.name

    def write(self, mpc_observation):
        """
        Writes a single entry in the Minor Planet Center's format.
        :param mpc_observation:
        """
        assert isinstance(mpc_observation, ObsRecord)
        try:
            key = mpc_observation.comment.frame.strip()
        except:
            key = mpc_observation.date.mjd

        try:
            # keep any 'discovery' flags that are set on observation if already in buffer.
            mpc_observation.discovery = self.buffer[key].discovery.is_discovery
        except Exception as ex:
            logging.debug(type(ex))
            logging.debug(str(ex))
            pass
        self.buffer[key] = mpc_observation

        if self.auto_flush:
            self.flush()

    def flush(self):
        for obs in self.get_chronological_buffered_observations():
            self._flush_observation(obs)

        self.filehandle.flush()

    def _flush_observation(self, obs):
        assert isinstance(obs, ObsRecord)
        if (self.auto_discovery and
                not obs.null_observation and
                not self._discovery_written):
            obs.discovery = True

        if obs.discovery:
            if self._discovery_written:
                obs.discovery.is_initial_discovery = False
            else:
                self._discovery_written = True

        if obs.date.jd not in self._written_mpc_observations:
            self._written_mpc_observations.append(obs.date.jd)
            line = self.formatter(obs)
            self.filehandle.write(line + "\n")

    def close(self):
        self.filehandle.close()

    def get_chronological_buffered_observations(self):
        jds = list(self.buffer.keys())
        jds.sort()
        sorted_obs = []
        for jd in jds:
            sorted_obs.append(self.buffer[jd])
        return sorted_obs


def make_tnodb_header(observations,
                      observatory_code=None,
                      observers=DEFAULT_OBSERVERS,
                      measurers=DEFAULT_MEASURERS,
                      telescope=DEFAULT_TELESCOPE,
                      astrometric_network=DEFAULT_ASTROMETRIC_NETWORK):
    """
    Write a header appropriate for a tnodb style of file.
    """
    observatory_code = observatory_code is None and observations[0].observatory_code or observatory_code

    odates = [obs.date for obs in observations]
    mindate = min(odates).iso.replace('-', '')[0:8]
    maxdate = max(odates).iso.replace('-', '')[0:8]

    header = "COD {}\n".format(observatory_code)

    if observers is not None:
        header += "OBS {}".format(observers[0])
        if len(observers) > 2:
            header += ", {}".format(", ".join(observers[1:-1]))
        if len(observers) > 1:
            header += " and {}".format(observers[-1])
        header += "\n"
    if measurers is not None:
        header += "MEA {}".format(measurers[0])
        if len(measurers) > 2:
            header += ", {}".format(", ".join(measurers[1:-1]))
        if len(measurers) > 1:
            header += " and {}".format(measurers[-1])
        header += "\n"
    header += "TEL {}\n".format(telescope)
    header += "NET {}\n".format(astrometric_network)
    header += mindate is not None and "{:s} {:s}\n".format('STD', mindate) or ""
    header += maxdate is not None and "{:s} {:s}\n".format('END', maxdate) or ""

    return header


class EphemerisReader(object):
    """
    A class to read in MPC files.

    Can be initialized with a filename, will then initialize the mpc_observations attribute to hold the observations.
    """

    def __init__(self, filename=None, replace_provisional=None, provisional_name=None):
        self.replace_provisional = replace_provisional
        self._provisional_name = provisional_name
        if filename is not None:
            self.filename = filename
            self.mpc_observations = self.read(filename)

    def read(self, filename):
        """
        Read  MPC records from filename:

        :param filename: filename of file like object.
        :rtype : [ObsRecord]
        """

        self.filename = filename
        # can be a file like objects,

        input_mpc_lines = None
        while True:
            try:
                if isinstance(filename, str):
                    filehandle = open(filename, "r")
                else:
                    filehandle = filename

                input_mpc_lines = filehandle.read().split('\n')
                filehandle.close()
                break
            except IOError as ioe:
                if ioe.errno == 16:
                    # Resource busy, so try again in a second.
                    import time
                    time.sleep(1)
                else:
                    raise ioe

        mpc_observations = []
        next_comment = None
        if input_mpc_lines is None:
            logging.warning("Failed to read any lines from file: {}".format(self.filename))
            return numpy.array([])

        for line in input_mpc_lines:
            line = line.rstrip()
            mpc_observation = ObsRecord.from_string(line)
            if isinstance(mpc_observation, OSSOSComment):
                next_comment = mpc_observation
                continue
            if isinstance(mpc_observation, ObserverLocation):
                mpc_observations[-1].location = mpc_observation
                continue
            if isinstance(mpc_observation, ObsRecord):
                if next_comment is not None:
                    mpc_observation.comment = next_comment
                    next_comment = None

                if self.replace_provisional is not None:  # then it has an OSSOS designation: set that in preference
                    mpc_observation.provisional_name = self.provisional_name
                if len(str(mpc_observation.provisional_name.strip())) == 0 and \
                        str(mpc_observation.minor_planet_number) == "":
                    mpc_observation.provisional_name = filename.split(".")[0]
                mpc_observations.append(mpc_observation)

        # No assurance that a .ast file is date-ordered: date-ordered is more expected behaviour
        mpc_observations.sort(key=lambda x: x.date)

        return numpy.array(mpc_observations)

    @property
    def provisional_name(self):
        """
        Determine the provisional name based on the file being accessed.
        :return: str
        """
        if self._provisional_name is not None:
            return self._provisional_name
        if isinstance(self.filename, str):
            self._provisional_name = self.filename.rstrip('.ast')
        elif hasattr(self.filename, 'name'):
            self._provisional_name = self.filename.name
        elif hasattr(self.filename, 'filename'):
            self._provisional_name = self.filename.filename
        elif hasattr(self.filename, '__class__'):
            self._provisional_name = str(self.filename.__class__)
        else:
            self._provisional_name = str(type(self.filename))
        self._provisional_name = os.path.basename(self._provisional_name)
        return self._provisional_name


class Index(object):
    """
    MOP/OSSOS name mapping index.
    """
    MAX_NAME_LENGTH = 10

    def __init__(self, idx_filename):
        self.names = {}
        self.index = {}
        if os.access(idx_filename, os.F_OK):
            with open(idx_filename, 'r') as idx_handle:
                for line in idx_handle.readlines():
                    master_name = line[0:Index.MAX_NAME_LENGTH]
                    master_name = master_name.strip()
                    self.names[master_name] = master_name
                    self.index[master_name] = [master_name]
                    for i in range(Index.MAX_NAME_LENGTH, len(line), Index.MAX_NAME_LENGTH):
                        this_name = line[i:i + Index.MAX_NAME_LENGTH].strip()
                        self.index[master_name].append(this_name)
                        self.names[this_name] = master_name

    def __format__(self, format_spec):
        return self.__str__()

    def __str__(self):
        result = ""
        for name in self.index:
            result += "{0:<{1}s}".format(name, Index.MAX_NAME_LENGTH)
            for alias in self.get_aliases(name):
                result += "{0:<{1}s}".format(alias, Index.MAX_NAME_LENGTH)
            result += "\n"
        return result

    def get_aliases(self, name):
        """
        get all names associated with a given name.
        :rtype : list
        :param name: object to get alias names of.
        """
        if name not in self.names:
            return name
        return self.index[self.names[name]]

    def is_same(self, name1, name2):
        """
        Do name1 and name2 refer to the same object?

        :param name1: name of object 1
        :param name2: name of object 2
        :return: Bool
        """
        return name2 in self.get_aliases(name1)


class MPCConverter(object):
    """
    Converts an MPC formatted file to a TNOdb one.
    :param mpc_file The input filename, of MPC lines.
    :param output   if required; else will use root of provided MPC file.

    batch_convert is factory method that will write out a series of input files given an input path.
    """

    def __init__(self, mpc_file, output=None):

        if output is None:
            output = mpc_file.rpartition('.')[0] + '.tnodb'

        self.mpc_file = mpc_file
        self.outfile = open(output, 'w')
        self.write_header = True

    def convert(self):
        with open(self.mpc_file, 'r') as infile:
            observations = []
            for line in infile.readlines():
                obs = ObsRecord().from_string(line)
                observations.append(obs)

            if self.write_header:
                self.outfile.write(make_tnodb_header(observations))
                self.write_header = False

            for obs in observations:
                self.outfile.write(obs.to_tnodb() + '\n')

    @classmethod
    def batch_convert(cls, path):
        for fn in os.listdir(path):
            if fn.endswith('.mpc') or fn.endswith('.track') or fn.endswith('.checkup') or fn.endswith('.nailing'):
                cls(path + fn).convert()


class CFEPSComment(OSSOSComment):
    """
    This holds the old-style comments that come for CFEPS style entries.
    """

    def __init__(self, frame, comment):

        if "measured inside confirm @" in comment:
            values = comment.split('@')[1].split()
            x = values[0]
            y = values[1]
        else:
            x = ""
            y = ""
        source_name = None
        mpc_note = " "
        super(CFEPSComment, self).__init__("L", frame, source_name, " ", mpc_note, x, y, comment=comment)

    @classmethod
    def from_string(cls, comment):
        """
        Build a comment from a CFEPS style comment string.
        """
        values = comment.split()
        if values[0] != "L" or len(values) < 2:
            raise ValueError("Not a CFEPS style comment: {}".format(comment))
        frame = values[1]
        comment = " ".join(values[2:])
        return cls(frame, comment)


class TNOdbComment(object):
    """
    This holds a TNOdb style comment line which contains flags.

    A TNOdb style comment consists of three space seperated fields that are used by the tnodb followed by a
    comment string that is either in the CFEPS or OSSOS format.
    """

    def __init__(self, index, date, flags, **kwargs):

        self.comment_object = None
        self.index = index
        self.date = date
        self.flags = flags

    @classmethod
    def from_string(cls, line):
        if len(line) < 56:
            raise ValueError("Line too short, not a valid TNOdb comment string: {}".format(line))
        index = line[0:14].strip()
        if not re.match(r'\d{8}_\S{3}_\S', index):
            raise ValueError("No valid flags, not a valid TNOdb comment string: {}".format(line))

        date = line[15:23].strip()
        flags = line[24:34].strip()
        comment = line[56:].strip()
        this = cls(index, date, flags)
        logging.debug("Got TNOdb parts: date: {} flags: {} comment: {}".format(date, flags, comment))
        comment_object = None
        # try build a comment object based on the TNOdb comment string
        if len(comment) > 0:
            for func in [OSSOSComment.from_string,
                         CFEPSComment.from_string,
                         MPCComment.from_string]:
                try:
                    comment_object = func(comment)
                except ValueError as verr:
                    logging.debug(verr)
                    continue
                break
        comment_object.date = date
        comment_object.flags = flags
        comment_object.index = index
        return comment_object


class RealOSSOSComment(OSSOSComment):
    @classmethod
    def from_string(cls, comment):
        if comment.strip()[0] != "O":
            comment = "O " + comment
        return super(RealOSSOSComment, cls).from_string(comment)

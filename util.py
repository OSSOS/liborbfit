"""OSSOS helper methods"""
from datetime import datetime
from logging import handlers
import logging
import numpy
import os
import re
import subprocess
import time
import vos
import sys

try:
    from astropy.time import erfa_time
except ImportError:
    from astropy.time import sofa_time as erfa_time
from astropy.time import TimeString
from astropy.time import Time


MATCH_TOLERANCE = 100.0


def task():
    return os.path.splitext(os.path.basename(sys.argv[0]))[0]


def exec_prog(args):
    """Run a subprocess, check for .OK and raise error if does not exist.

    args:  list of arguments, for value is the command to execute.
    """

    program_name = args[0]
    logging.info(" ".join(args))
    output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    if not os.access(program_name+".OK", os.F_OK):
        logging.error("No {}.OK file?".format(program_name))
        raise subprocess.CalledProcessError(-1, ' '.join(args), output)
    os.unlink(program_name+".OK")
    if os.access(program_name+".FAILED", os.F_OK):
        os.unlink(program_name+".FAILED")
    return output


class VOFileHandler(handlers.BufferingHandler):
    """
    A handler class that writes formatted logging records to VOSpace files.
    """
    def __init__(self, filename, vos_client=None):
        self.filename = filename
        self._client = vos_client
        self.stream = None
        super(VOFileHandler, self).__init__(1024*1024)

    @property
    def client(self):
        """
        Send back the client we were sent, or construct a default one.
        """
        if self._client is not None:
            return self._client
        self._client = vos.Client()
        return self._client

    def close(self):
        """
        Closes the stream.
        """
        self.flush()
        try:
            if self.stream is not None:
                self.stream.close()
                self.stream = None
        except:
            pass

    def flush(self):
        if self.stream is None:
            self.stream = self._open()
        for record in self.buffer:
            self.stream.write("{}\n".format(self.format(record)))
        self.buffer = []

    def _open(self):
        """
        Open the current base file with the (original) mode
        Return the resulting stream.
        """
        try:
            if self.stream.closed:
                return self.stream
        except AttributeError:
            pass

        buff = None
        try:
            buff = self.client.open(self.filename, view='data').read()
        except:
            pass

        stream = self.client.open(self.filename, os.O_WRONLY)
        if buff is not None:
            stream.write(buff)

        return stream


def get_pixel_bounds_from_datasec_keyword(datasec):
    """
    Return the x/y pixel boundaries of the data section.
    :param datasec: str e.g. '[33:2080,1:4612]'
    :return: ((xmin,xmax),(ymin,ymax))
    """
    datasec = re.findall(r'(\d+)', datasec)
    x1 = min(int(datasec[0]), int(datasec[1]))
    x2 = max(int(datasec[0]), int(datasec[1]))
    y1 = min(int(datasec[2]), int(datasec[3]))
    y2 = max(int(datasec[2]), int(datasec[3]))

    return (x1, x2), (y1, y2)


def match_lists(pos1, pos2, tolerance=MATCH_TOLERANCE):
    """
    Given two sets of x/y positions match the lists, uniquely.

    :rtype : numpy.ma, numpy.ma
    :param pos1: list of x/y positions.
    :param pos2: list of x/y positions.
    :param tolerance: float distance, in pixels, to consider a match

    Algorithm:
        - Find all the members of pos2 that are within tolerance of pos1[idx1].
                These pos2 members are match_group_1
        - Find all the members of pos1 that are within tolerance of match_group_1[idx2].
                These pos1 members are match_group_2
        - If pos1[idx] is in match_group_2 then pos1[idx] is a match of object at match_group_1[idx2]

    """

    assert isinstance(pos1, numpy.ndarray)
    assert isinstance(pos2, numpy.ndarray)

    # build some arrays to hold the index of things that matched between lists.
    npts1 = len(pos1[:, 0])
    pos1_idx_array = numpy.arange(npts1, dtype=numpy.int16)
    npts2 = len(pos2[:, 0])
    pos2_idx_array = numpy.arange(npts2, dtype=numpy.int16)

    # this is the array of final matched index, -1 indicates no match found.
    match1 = numpy.ma.zeros(npts1, dtype=numpy.int16)
    match1.mask = True

    # this is the array of matches in pos2, -1 indicates no match found.
    match2 = numpy.ma.zeros(npts2, dtype=numpy.int16)
    match2.mask = True

    for idx1 in range(npts1):

        # compute the distance source idx1 to each member of pos2
        sep = numpy.sqrt((pos2[:, 0] - pos1[idx1, 0]) ** 2 + (pos2[:, 1] - pos1[idx1, 1]) ** 2)

        # considered a match if sep is below tolerance and is the closest match available.
        match_condition = numpy.all((sep <= tolerance, sep == sep.min()), axis=0)

        # match_group_1 is list of the indexes of pos2 entries that qualified as possible matches to pos1[idx1]
        match_group_1 = pos2_idx_array[match_condition]

        # For each of those pos2 objects that could be a match to pos1[idx] find the best match in all of pos1
        for idx2 in match_group_1:
            # compute the distance from this pos2 object that is a possible match to pos1[idx1] to all members of pos1
            sep = numpy.sqrt((pos1[:, 0] - pos2[idx2, 0]) ** 2 + (pos1[:, 1] - pos2[idx2, 1]) ** 2)

            # considered a match if sep is below tolerance and is the closest match available.
            match_condition = numpy.all((sep <= tolerance, sep == sep.min()), axis=0)
            match_group_2 = pos1_idx_array[match_condition]

            # Are any of the pos1 members that were matches to the matched pos2 member the pos1[idx] entry?
            if idx1 in match_group_2:
                match1[idx1] = idx2
                match2[idx2] = idx1
                # this BREAK is in here since once we have a match we're done.
                break

    return match1, match2


class TimeMPC(TimeString):
    """
    Override the TimeString class to convert from MPC format string to astropy.time.Time object.

    usage:

    from astropy.time.core import Time
    Time.FORMATS[TimeMPC.name] = TimeMPC

    t = Time('2000 01 01.00001', format='mpc', scale='utc')

    str(t) == '2000 01 01.000001'
    """

    name = 'mpc'
    subfmts = (('mpc', '%Y %m %d', "{year:4d} {mon:02d} {day:02d}.{fracday:s}"),)

    def __init__(self, val1, val2, scale, precision,
                 in_subfmt, out_subfmt, from_jd=False):
        super(TimeMPC, self).__init__(val1=val1,
                                      val2=val2,
                                      scale=scale,
                                      precision=precision,
                                      in_subfmt=in_subfmt,
                                      out_subfmt=out_subfmt,
                                      from_jd=from_jd)
        self.precision = precision

    # ## need our own 'set_jds' function as the MPC Time string is not typical
    def set_jds(self, val1, val2):
        """

        Parse the time strings contained in val1 and set jd1, jd2

        :param val1: array of strings to parse into JD format
        :param val2: not used for string conversions but passed regardless
        """

        # This routine is based on atropy.time.core.TimeString class.

        iterator = numpy.nditer([val1, None, None, None, None, None, None],
                                op_dtypes=[val1.dtype] + 5 * [numpy.intc] + [numpy.double])
        subfmts = self.subfmts
        for val, iy, im, iday, ihr, imin, dsec in iterator:
            time_str = val.item()

            # Assume that anything following "." on the right side is a
            # floating fraction of a day.
            try:
                idot = time_str.rindex('.')
            except:
                fracday = 0.0
            else:
                time_str, fracday = time_str[:idot], time_str[idot:]
                fracday = float(fracday)

            for _, strptime_fmt, _ in subfmts:
                try:
                    tm = time.strptime(time_str, strptime_fmt)
                except ValueError:
                    pass
                else:
                    iy[...] = tm.tm_year
                    im[...] = tm.tm_mon
                    iday[...] = tm.tm_mday
                    hrprt = tm.tm_hour + int(24 * fracday)
                    ihr[...] = hrprt
                    mnprt = tm.tm_min + int(60 * (24 * fracday - hrprt))
                    imin[...] = mnprt
                    dsec[...] = tm.tm_sec + 60 * (60 * (24 * fracday - hrprt) - mnprt)

                    break
            else:
                raise ValueError("Time {0} does not match {1} format".format(time_str, self.name))

        self.jd1, self.jd2 = erfa_time.dtf_jd(
            self.scale.upper().encode('utf8'), *iterator.operands[1:])

        return

    def str_kwargs(self):
        """
        Generator that yields a dict of values corresponding to the
        calendar date and time for the internal JD values.
        """
        iys, ims, ids, ihmsfs = erfa_time.jd_dtf(self.scale.upper()
                                                 .encode('utf8'),
                                                 6,
                                                 self.jd1, self.jd2)

        # Get the str_fmt element of the first allowed output subformat
        _, _, str_fmt = self._select_subfmts(self.out_subfmt)[0]

        yday = None
        has_yday = '{yday:' in str_fmt or False

        ihrs = ihmsfs[..., 0]
        imins = ihmsfs[..., 1]
        isecs = ihmsfs[..., 2]
        ifracs = ihmsfs[..., 3]
        for iy, im, iday, ihr, imin, isec, ifracsec in numpy.nditer(
                [iys, ims, ids, ihrs, imins, isecs, ifracs]):
            if has_yday:
                yday = datetime(iy, im, iday).timetuple().tm_yday

            fracday = (((((ifracsec / 1000000.0 + isec) / 60.0 + imin) / 60.0) + ihr) / 24.0) * (10 ** 6)
            fracday = '{0:06g}'.format(fracday)[0:self.precision]

            yield {'year': int(iy), 'mon': int(im), 'day': int(iday),
                   'hour': int(ihr), 'min': int(imin), 'sec': int(isec),
                   'fracsec': int(ifracsec), 'yday': yday, 'fracday': fracday}

Time.FORMATS['mpc'] = TimeMPC

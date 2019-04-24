import logging
from datetime import datetime
import numpy
import re
import time
import six
from astropy._erfa import d2dtf, dtf2d
from astropy.time import TimeString


class TimeMPC(TimeString):
    """
    Override the TimeString class to convert from MPC format string to astropy.time.Time object.

    usage:

    from time_mpc import TimeMPC
    from astropy.time.core import Time

    t = Time('2000 01 01.00001', format='mpc', scale='utc')

    str(t) == '2000 01 01.000001'
    """
    name = 'mpc'
    subfmts = (('mpc', '%Y %m %d', "{year:4d} {mon:02d} {day:02d}.{fracday:s}"),)

    def __init__(self, val1, val2, scale, precision=6,
                 in_subfmt=None, out_subfmt=None, from_jd=False):
        super(TimeMPC, self).__init__(val1=val1,
                                      val2=val2,
                                      scale=scale,
                                      precision=precision,
                                      in_subfmt=in_subfmt,
                                      out_subfmt=out_subfmt,
                                      from_jd=from_jd)
        self.precision = precision

    def set_jds(self, val1, val2):
        """Parse the time strings contained in val1 and set jd1, jd2"""
        # Select subformats based on current self.in_subfmt
        subfmts = self._select_subfmts(self.in_subfmt)

        iterator = numpy.nditer([val1, None, None, None, None, None, None],
                                op_dtypes=[val1.dtype] + 5 * [numpy.intc] + [numpy.double])

        for val, iy, im, i_day, ihr, i_min, d_sec in iterator:
            iy[...], im[...], i_day[...], ihr[...], i_min[...], d_sec[...] = (
                self.parse_string(val.item(), subfmts))

        self.jd1, self.jd2 = dtf2d(
            self.scale.upper().encode('utf8'), *iterator.operands[1:])

    def parse_string(self, timestr, subfmts):
        """
        Read time from a single string, using a set of possible formats.

        :param timestr:  A string representing the time. eg:  '2001 01 01.00001'
        :param subfmts: which format is the string in, eg: 'mpc'
        """
        # Datetime components required for conversion to JD by ERFA, along
        # with the default values.
        components = ('year', 'mon', 'mday')
        defaults = (None, 1, 1, 0)
        # Assume that anything following "." on the right side is a
        # floating fraction of a second.
        try:
            idot = timestr.rindex('.')
        except:
            fracday = 0.0
        else:
            timestr, fracday = timestr[:idot], timestr[idot:]
            fracday = float(fracday)

        for _, strptime_fmt_or_regex, _ in subfmts:
            if isinstance(strptime_fmt_or_regex, six.string_types):
                try:
                    tm = time.strptime(timestr, strptime_fmt_or_regex)
                except ValueError as ex:
                    logging.debug(str(ex))
                    continue
                else:
                    vals = [getattr(tm, 'tm_' + component)
                            for component in components]
                    vals.append(tm.tm_hour + int(24 * fracday))
                    vals.append(tm.tm_min + int(60 * (24 * fracday - vals[-1])))
                    vals.append(tm.tm_sec + 60 * (60 * (24 * fracday - vals[-2]) - vals[-1]))
            else:
                tm = re.match(strptime_fmt_or_regex, timestr)
                if tm is None:
                    continue
                tm = tm.groupdict()
                vals = [int(tm.get(component, default)) for component, default
                        in six.moves.zip(components, defaults)]

                hrprt = int(24 * fracday)
                vals.append(hrprt)
                mnprt = int(60 * (24 * fracday - hrprt))
                vals.append(mnprt)
                scprt = 60 * (60 * (24 * fracday - hrprt) - mnprt)
                vals.append(scprt)
            return vals
        else:
            raise ValueError('Time {0} does not match {1} format'
                             .format(timestr, self.name))

    def str_kwargs(self):
        """
        Generator that yields a dict of values corresponding to the
        calendar date and time for the internal JD values.
        """
        scale = self.scale.upper().encode('ascii'),
        iys, ims, ids, ihmsfs = d2dtf(scale, self.precision,
                                      self.jd1, self.jd2_filled)

        # Get the str_fmt element of the first allowed output subformat
        _, _, str_fmt = self._select_subfmts(self.out_subfmt)[0]

        yday = None
        has_yday = True if '{yday:' in str_fmt else False

        ihrs = ihmsfs['h']
        imins = ihmsfs['m']
        isecs = ihmsfs['s']
        ifracs = ihmsfs['f']
        for iy, im, id, ihr, imin, isec, ifracsec in numpy.nditer(
                [iys, ims, ids, ihrs, imins, isecs, ifracs]):
            if has_yday:
                yday = datetime(iy, im, id).timetuple().tm_yday

            fracday = (((((ifracsec / 1000000.0 + isec) / 60.0 + imin) / 60.0) + ihr) / 24.0) * (10 ** 6)
            fracday = '{0:06g}'.format(fracday)[0:self.precision]

            yield {'year': int(iy), 'mon': int(im), 'day': int(id),
                   'hour': int(ihr), 'min': int(imin), 'sec': int(isec),
                   'fracsec': int(ifracsec), 'yday': yday, 'fracday': fracday}

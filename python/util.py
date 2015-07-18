"""OSSOS helper methods"""
from datetime import datetime
import numpy
import time

try:
    from astropy.time import erfa_time
except ImportError:
    from astropy.time import sofa_time as erfa_time
from astropy.time import TimeString
from astropy.time import Time


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

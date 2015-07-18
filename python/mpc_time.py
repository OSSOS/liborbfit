import itertools

__author__ = 'jjk'

from datetime import datetime
import time
import numpy

from astropy.time import erfa_time as sofa_time
from astropy.time import TimeString
from astropy.time import Time

class TimeMPC(TimeString):
    """
    Override the TimeString class to convert from MPC format string to astropy.time.Time object.

    usage:

    from astropy.time import Time
    Time.FORMATS[TimeMPC.name] = TimeMPC

    t = Time('2000 01 01.00001', format='mpc', scale='utc')

    str(t) == '2000 01 01.000001'
    """

    name = 'mpc'
    subfmts = (('mpc', '%Y %m %d', "{year:4d} {mon:02d} {day:02d}.{fracday:s}"),)

    ### need our own 'set_jds' function as the MPC Time string is not typical
    def set_jds(self, val1, val2):
        """

        Parse the time strings contained in val1 and set jd1, jd2

        :param val1: array of strings to parse into JD format
        :param val2: not used for string conversions but passed regardless
        """
        n_times = len(val1)  # val1,2 already checked to have same len
        iy = numpy.empty(n_times, dtype=numpy.intc)
        im = numpy.empty(n_times, dtype=numpy.intc)
        iday = numpy.empty(n_times, dtype=numpy.intc)
        ihr = numpy.empty(n_times, dtype=numpy.intc)
        imin = numpy.empty(n_times, dtype=numpy.intc)
        dsec = numpy.empty(n_times, dtype=numpy.double)

        # Select subformats based on current self.in_subfmt
        subfmts = self._select_subfmts(self.in_subfmt)

        for i, timestr in enumerate(val1):
            # Assume that anything following "." on the right side is a
            # floating fraction of a day.
            try:
                idot = timestr.rindex('.')
            except:
                fracday = 0.0
            else:
                timestr, fracday = timestr[:idot], timestr[idot:]
                fracday = float(fracday)

            for _, strptime_fmt, _ in subfmts:
                try:
                    tm = time.strptime(timestr, strptime_fmt)
                except ValueError:
                    pass
                else:
                    iy[i] = tm.tm_year
                    im[i] = tm.tm_mon
                    iday[i] = tm.tm_mday
                    ihr[i] = tm.tm_hour + int(24 * fracday)
                    imin[i] = tm.tm_min + int(60 * (24 * fracday - ihr[i]))
                    dsec[i] = tm.tm_sec + 60 * (60 * (24 * fracday - ihr[i]) - imin[i])
                    break
            else:
                raise ValueError('Time {0} does not match {1} format'
                .format(timestr, self.name))

        self.jd1, self.jd2 = sofa_time.dtf_jd(self.scale.upper().encode('utf8'),
                                              iy, im, iday, ihr, imin, dsec)
        return

    def str_kwargs(self):
        """                                                                                                                                           
        Generator that yields a dict of values corresponding to the                                                                                   
        calendar date and time for the internal JD values.

        Here we provide the additional 'fracday' element needed by 'mpc' format
        """
        iys, ims, ids, ihmsfs = sofa_time.jd_dtf(self.scale.upper()
                                                 .encode('utf8'),
                                                 6,
                                                 self.jd1, self.jd2)

        # Get the str_fmt element of the first allowed output subformat                                                                               
        _, _, str_fmt = self._select_subfmts(self.out_subfmt)[0]

        yday = None
        has_yday = '{yday:' in str_fmt or False

        for iy, im, iday, ihmsf in itertools.izip(iys, ims, ids, ihmsfs):
            ihr, imin, isec, ifracsec = ihmsf
            if has_yday:
                yday = datetime(iy, im, iday).timetuple().tm_yday

            # MPC uses day fraction as time part of datetime
            fracday = (((((ifracsec / 1000000.0 + isec) / 60.0 + imin) / 60.0) + ihr) / 24.0) * (10 ** 6)
            fracday = '{0:06g}'.format(fracday)[0:self.precision]
            #format_str = '{{0:g}}'.format(self.precision)
            #fracday = int(format_str.format(fracday))
            yield dict(year=int(iy), mon=int(im), day=int(iday), hour=int(ihr), min=int(imin), sec=int(isec),
                       fracsec=int(ifracsec), yday=yday, fracday=fracday)


Time.FORMATS[TimeMPC.name] = TimeMPC


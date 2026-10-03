"""MPC calendar-date time format for :class:`astropy.time.Time`.

Registers format ``mpc`` so strings like ``YYYY MM DD.ffffff`` (where the
fraction is a fraction of a day) can initialize and format ``Time`` objects.

Subclassing :class:`~astropy.time.TimeString` is enough for Astropy 6+: the
base class owns JD conversion, masking, and output. This module only teaches
it that digits after the decimal are a **day** fraction (not a second
fraction, which is the ``TimeString`` default).
"""

import re

from astropy.time import TimeString

# Importing this module registers the format via TimeString's metaclass.
__all__ = ["TimeMPC"]


class TimeMPC(TimeString):
    """
    Minor Planet Center calendar date: ``YYYY MM DD.ffffff``.

    The fractional part is a fraction of a **day** (MPC convention), not a
    fraction of a second.

    Examples
    --------
    >>> from astropy.time import Time  # doctest: +SKIP
    >>> import mp_ephem.time_mpc  # registers format='mpc'  # doctest: +SKIP
    >>> t = Time('2000 01 01.00001', format='mpc', scale='utc', precision=5)
    >>> t.iso
    '2000-01-01 00:00:00.86400'
    >>> t.mpc
    '2000 01 01.00001'
    """

    name = "mpc"
    subfmts = (
        (
            "mpc",
            re.compile(
                r"^\s*(?P<year>\d{4})\s+(?P<mon>\d{1,2})\s+(?P<mday>\d{1,2})"
                r"(?:\.(?P<fracday>\d+))?\s*$"
            ),
            "{year:4d} {mon:02d} {day:02d}.{fracday:s}",
        ),
    )

    def parse_string(self, timestr, subfmts):
        """Parse one MPC date string into ``(year, mon, mday, hour, min, sec)``."""
        for _, parser, _ in subfmts:
            match = parser.match(timestr)
            if match is None:
                continue
            groups = match.groupdict()
            year = int(groups["year"])
            mon = int(groups["mon"])
            mday = int(groups["mday"])
            frac_digits = groups["fracday"] or ""
            fracday = (
                int(frac_digits) / (10 ** len(frac_digits)) if frac_digits else 0.0
            )
            seconds = fracday * 86400.0
            hour = int(seconds // 3600)
            minute = int((seconds % 3600) // 60)
            sec = seconds - hour * 3600 - minute * 60
            return [year, mon, mday, hour, minute, sec]
        raise ValueError(f"Time {timestr} does not match {self.name} format")

    def str_kwargs(self):
        """Yield format kwargs, adding MPC ``fracday`` from h/m/s/fracsec."""
        for kwargs in super().str_kwargs():
            total_sec = (
                (kwargs["hour"] * 60 + kwargs["min"]) * 60
                + kwargs["sec"]
                + kwargs["fracsec"] / (10 ** self.precision)
            )
            if self.precision > 0:
                frac_int = int(round(total_sec / 86400.0 * (10 ** self.precision)))
                modulus = 10 ** self.precision
                if frac_int >= modulus:
                    frac_int -= modulus
                kwargs["fracday"] = f"{frac_int:0{self.precision}d}"
            else:
                kwargs["fracday"] = ""
            yield kwargs

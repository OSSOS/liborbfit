#!/usr/bin/env python

from mp_ephem.ephem_target import EphemTarget
from astropy.time import Time
from astropy import units
from mp_ephem import horizons
import ephem
from mp_ephem import bk_orbit
import math
from copy import deepcopy
import argparse
import sys
import os


def build_ephem_files(target_name, start_time, stop_time, step_size=None, observatory=None,
                      ephem_format=None, runid=None):

    if observatory is None:
        _cfht = ephem.Observer()
        _cfht.lat = 0.344
        _cfht.lon = -2.707
        _cfht.elevation = 4100
        _cfht.date = '2018/03/28 20:00:00'
        observatory = _cfht
    if step_size is None:
        step_size = 30*units.minute
    start_time = Time(start_time)
    stop_time = Time(stop_time)
    sun = ephem.Sun()
    fb = ephem.FixedBody()


    if os.access(target_name, os.F_OK):
        body = bk_orbit.BKOrbit(None, ast_filename=target_name)
        target_name = body.name
    else:
        body = horizons.Body(target_name.replace("_", " "), start_time=start_time, stop_time=stop_time,
                         step_size=step_size, center='568')
    current_time = start_time
    body.predict(current_time)

    et = EphemTarget(target_name.replace(" ", "_"), ephem_format=ephem_format, runid=runid)

    while current_time < stop_time:
        target_up = False
        observatory.date = current_time.iso.replace('-', '/')
        observatory.horizon = math.radians(-7)
        sun.compute(observatory)
        sun_rise = Time(str(sun.rise_time).replace('/', '-'))
        sun_set = Time(str(sun.set_time).replace('/', '-'))

        if current_time < sun_set or current_time > sun_rise:
            sun_down = False
        else:
            sun_down = True

        body.predict(current_time)
        fb._ra = body.coordinate.ra.radian
        fb._dec = body.coordinate.dec.radian

        observatory.horizon = math.radians(40)
        fb.compute(observatory)
        fb_rise_time = Time(str(fb.rise_time).replace('/', '-'))
        fb_set_time = Time(str(fb.set_time).replace('/', '-'))

        if fb_rise_time < current_time < fb_set_time:
            target_up = True
        if fb_rise_time < current_time and fb_set_time < current_time:
            target_up = False
        if fb_rise_time > current_time and fb_set_time > current_time:
            target_up = True
        if fb_rise_time > current_time > fb_set_time:
            target_up = False

        if target_up and sun_down:
            coordinate = deepcopy(body.coordinate)
            coordinate.mag = body.mag
            coordinate.obstime = current_time
            et.append(coordinate)
        current_time += step_size

    et.save()
    print coordinate.mag


def minor_planet_ephem(target_names, start_time, stop_time, step_size=None, observatory=None, ephem_format=None, runid=None):
    """
    Given a list of targets build an ephemeris file to load to CFHT
    This routine will only put out lines for when the target is up.

    :param target_names:
    :param start_time:
    :param stop_time:
    :param step_size:
    :param observatory:
    :return:
    """

    start_time = Time(start_time)
    stop_time = Time(stop_time)

    if step_size is None:
        step_size = 30 * units.minute

    for target_name in target_names:
        build_ephem_files(target_name, start_time, stop_time, step_size=step_size, observatory=observatory,
                          ephem_format=ephem_format, runid=runid)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('start_time', help="Date at start of dark run.")
    parser.add_argument('end_time', help="Date at end of dark run.")
    parser.add_argument('target_names', nargs="+", help="Names of targets to build ephemeris files for.")
    parser.add_argument('--runid', default='17AC99')
    parser.add_argument('--ephem-format', default='CFHT ET')
    parser.add_argument('--step-size', help="size of time step for ephemeris.", default=300 * units.minute)
    parser.add_argument('--observatory', default=None)

    args = parser.parse_args()
    if not isinstance(args.step_size, units.Quantity): 
       args.step_size *= units.minute
    minor_planet_ephem(args.target_names, args.start_time, args.end_time, args.step_size, args.observatory, args.ephem_format,
         args.runid)

if __name__ == '__main__':
    sys.exit(main())

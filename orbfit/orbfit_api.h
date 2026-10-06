/* Entry points of the orbfit shared library, as called from mp_ephem/bk_orbit.py
 * through ctypes.  Each returns a pointer to a static result array that is
 * overwritten by the next call; the library is not reentrant.
 *
 * On failure the results are NaN and orbfit_last_error() describes the
 * problem; it is the empty string after a successful call.
 *
 * The library is built with hidden symbol visibility; only the functions
 * marked ORBFIT_API here are exported. */
#ifndef ORBFIT_API_H
#define ORBFIT_API_H

#if defined(__GNUC__)
#define ORBFIT_API __attribute__((visibility("default")))
#else
#define ORBFIT_API
#endif

/* Fit the observations in mpc_filename and write the a/b/g orbit to
 * abg_filename.  Returns {distance, distance uncertainty} in AU. */
ORBFIT_API double *fitradec(char *mpc_filename, char *abg_filename);

/* RA, Dec (deg), error ellipse a, b (arcsec), PA (deg), distance (AU),
 * ecliptic lon, lat (deg) at UTC Julian Date jdate from site obscode. */
ORBFIT_API double *predict(char *abg_file, double jdate, int obscode);

/* Position {x, y, z} (AU) of the target at UTC Julian Date jdate. */
ORBFIT_API double *predict_helio(char *abg_file, double jdate, int obscode);

/* a, e, i, Node, peri, T, their uncertainties, epoch, distance and distance
 * uncertainty for the orbit in abg_file. */
ORBFIT_API double *abg_to_aei(char *abg_file);

/* Message from the last failed call, or "" */
ORBFIT_API const char *orbfit_last_error(void);

/* Override the ephemeris and observatory files (otherwise taken from the
 * ORBIT_EPHEMERIS and ORBIT_OBSERVATORIES environment variables) and the
 * default astrometric uncertainty (arcsec) of MPC-format observations. */
ORBFIT_API void set_ephem_file(char *fname);
ORBFIT_API void set_observatory_file(char *fname);
ORBFIT_API void set_mpc_dtheta(double d);

/* TT Julian Date for a UTC Julian Date */
ORBFIT_API double utc_jd_to_tt(double jd_utc);

#endif

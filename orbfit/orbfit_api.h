/* Entry points of the orbfit shared library, as called from mp_ephem/bk_orbit.py
 * through ctypes.  Each returns a pointer to a static result array that is
 * overwritten by the next call; the library is not reentrant. */
#ifndef ORBFIT_API_H
#define ORBFIT_API_H

/* Fit the observations in mpc_filename and write the a/b/g orbit to
 * abg_filename.  Returns {distance, distance uncertainty} in AU. */
double *fitradec(char *mpc_filename, char *abg_filename);

/* RA, Dec (deg), error ellipse a, b (arcsec), PA (deg), distance (AU),
 * ecliptic lon, lat (deg) at UTC Julian Date jdate from site obscode. */
double *predict(char *abg_file, double jdate, int obscode);

/* Position {x, y, z} (AU) of the target at UTC Julian Date jdate. */
double *predict_helio(char *abg_file, double jdate, int obscode);

/* a, e, i, Node, peri, T, their uncertainties, epoch, distance and distance
 * uncertainty for the orbit in abg_file. */
double *abg_to_aei(char *abg_file);

#endif

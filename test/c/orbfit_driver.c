/* Exercise the orbfit library entry points without Python, so the C code can
 * be run under AddressSanitizer / LeakSanitizer / UBSan.
 *
 *   orbfit_driver <observations> <abg output> [repeat]
 *
 * ORBIT_EPHEMERIS and ORBIT_OBSERVATORIES must point at the data files.
 * Exits non-zero if a returned position or element is not finite, or if a
 * call that should fail does not report an error.
 * Uncertainties may be NaN when the arc leaves the orbit undetermined. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "orbfit_api.h"

/* required[] lists the indices of v that must be finite, ending with -1 */
static int
check(const char *what, const double *v, int n, const int *required, int show)
{
  int i, bad = 0;
  if (show) {
    printf("%-14s", what);
    for (i = 0; i < n; i++) printf(" %.10g", v[i]);
    printf("\n");
  }
  for (i = 0; required[i] >= 0; i++)
    if (!isfinite(v[required[i]])) bad = 1;
  if (bad) fprintf(stderr, "non-finite value returned by %s\n", what);
  if (orbfit_last_error()[0]) {
    fprintf(stderr, "%s reported an error: %s\n", what, orbfit_last_error());
    bad = 1;
  }
  return bad;
}

/* A failed call must return NaN and leave a message */
static int
check_failed(const char *what, const double *v)
{
  if (!isnan(v[0]) || orbfit_last_error()[0] == 0) {
    fprintf(stderr, "%s did not report a failure\n", what);
    return 1;
  }
  return 0;
}

int
main(int argc, char *argv[])
{
  static const int fit_req[] = {0, -1};
  static const int aei_req[] = {0, 1, 2, 3, 4, 5, 12, 13, -1};
  static const int predict_req[] = {0, 1, 5, 6, 7, -1};
  static const int helio_req[] = {0, 1, 2, -1};
  double *r, jd;
  int i, repeat, bad = 0;

  if (argc < 3) {
    fprintf(stderr, "usage: %s <observations> <abg output> [repeat]\n", argv[0]);
    return 2;
  }
  repeat = argc > 3 ? atoi(argv[3]) : 1;

  for (i = 0; i < repeat; i++) {
    r = fitradec(argv[1], argv[2]);
    bad |= check("fitradec", r, 2, fit_req, i == 0);
    r = abg_to_aei(argv[2]);
    bad |= check("abg_to_aei", r, 15, aei_req, i == 0);
    jd = r[12];
    r = predict(argv[2], jd + 365.25, 568);
    bad |= check("predict", r, 8, predict_req, i == 0);
    r = predict_helio(argv[2], jd + 365.25, 568);
    bad |= check("predict_helio", r, 3, helio_req, i == 0);

    /* Failures part way through a calculation must not exit or leak, and
     * the next call must work. */
    bad |= check_failed("predict out of ephemeris range",
			predict(argv[2], 1.0e7, 568));
    bad |= check_failed("predict_helio out of ephemeris range",
			predict_helio(argv[2], 1.0e7, 568));
    bad |= check_failed("abg_to_aei missing file",
			abg_to_aei("/nonexistent/orbfit.abg"));
    bad |= check_failed("fitradec missing file",
			fitradec("/nonexistent/orbfit.mpc", argv[2]));
    r = predict(argv[2], jd + 365.25, 568);
    bad |= check("predict after failure", r, 8, predict_req, 0);
  }
  return bad;
}

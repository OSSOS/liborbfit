/* Exercise the orbfit library entry points without Python, so the C code can
 * be run under AddressSanitizer / LeakSanitizer / UBSan.
 *
 *   orbfit_driver <observations> <abg output> [repeat]
 *
 * ORBIT_EPHEMERIS and ORBIT_OBSERVATORIES must point at the data files.
 * Exits non-zero if any returned value is not finite. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "orbfit_api.h"

static int
check(const char *what, const double *v, int n, int show)
{
  int i, bad = 0;
  if (show) printf("%-14s", what);
  for (i = 0; i < n; i++) {
    if (!isfinite(v[i])) bad = 1;
    if (show) printf(" %.10g", v[i]);
  }
  if (show) printf("\n");
  if (bad) fprintf(stderr, "non-finite value returned by %s\n", what);
  return bad;
}

int
main(int argc, char *argv[])
{
  double *r, jd;
  int i, repeat, bad = 0;

  if (argc < 3) {
    fprintf(stderr, "usage: %s <observations> <abg output> [repeat]\n", argv[0]);
    return 2;
  }
  repeat = argc > 3 ? atoi(argv[3]) : 1;

  for (i = 0; i < repeat; i++) {
    r = fitradec(argv[1], argv[2]);
    bad |= check("fitradec", r, 2, i == 0);
    r = abg_to_aei(argv[2]);
    bad |= check("abg_to_aei", r, 15, i == 0);
    jd = r[12];
    r = predict(argv[2], jd + 365.25, 568);
    bad |= check("predict", r, 8, i == 0);
    r = predict_helio(argv[2], jd + 365.25, 568);
    bad |= check("predict_helio", r, 3, i == 0);
  }
  return bad;
}

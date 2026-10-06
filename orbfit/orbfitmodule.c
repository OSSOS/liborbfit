#include "orbfit.h"
#include "orbfit_api.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* Each entry point installs a guard so that orbfit_fail() anywhere below
 * returns here: the results are left as NaN, every Numerical Recipes
 * allocation still live is released and the integrator state is dropped. */
#define ORBFIT_GUARD(guard, result, n)		\
  do {						\
    fill_nan(result, n);			\
    if (setjmp(guard)) {			\
      nr_free_all();				\
      kbo3d_reset();				\
      fill_nan(result, n);			\
      orbfit_end();				\
      return result;				\
    }						\
    orbfit_begin(&guard);			\
  } while (0)

static void
fill_nan(double *result, int n)
{
  int i;
  for (i = 0; i < n; i++) result[i] = NAN;
}

static void
write_abg(const char *abg_filename, PBASIS *p, double **covar)
{
  FILE *abg_file;

  if ((abg_file = fopen(abg_filename,"w")) == NULL)
    orbfit_fail("Error opening a/b/g output file %s", abg_filename);

  fprintf(abg_file, "# Exact a, adot, b, bdot, g, gdot:\n");
  fprintf(abg_file, "%.17g %.17g %.17g %.17g %.17g %.17g\n",p->a,p->adot,p->b,
  	p->bdot, p->g, p->gdot);

  fprintf(abg_file, "# Covariance matrix: \n");
  print_matrix(abg_file, covar, 6, 6);

  /* Print out information on the coordinate system */
  fprintf(abg_file, "#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  fprintf(abg_file, "%.17g %.17g %.17g %.17g %.17g %.17g\n",
	 lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);

  if (fclose(abg_file) != 0)
    orbfit_fail("Error writing a/b/g output file %s", abg_filename);
}

double *fitradec(char *mpc_filename, char *abg_filename)
{
  static double result[2];
  jmp_buf guard;
  OBSERVATION *obsarray;
  int     nobs;
  PBASIS p;
  double d, dd;
  double **covar;
  double chisq;
  int dof;

  ORBFIT_GUARD(guard, result, 2);

  obsarray = nr_malloc(MAXOBS*sizeof(OBSERVATION));
  covar = dmatrix(1,6,1,6);

  if (read_radec(obsarray, mpc_filename, &nobs)) orbfit_rethrow();

  /* Call subroutine to do the actual fitting: */
  if (fit_observations(obsarray, nobs, &p, covar, &chisq, &dof, NULL) < 0)
    orbfit_rethrow();

  write_abg(abg_filename, &p, covar);

  /* Barycentric distance and its uncertainty */
  d = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
  dd = d*d*sqrt(covar[5][5]);
  result[0] = d;
  result[1] = dd;

  free_dmatrix(covar,1,6,1,6);
  nr_free(obsarray);
  orbfit_end();
  return result;
}

/* predict_helio.c - Read a file containing the a/b/g orbit fit and spit out
 * the predicted barycentric position
 * 9/20/19 jjk
 */

double *predict_helio(char *abg_file, double jdate, int obscode) {

  static double result[3];
  jmp_buf guard;
  PBASIS p;
  OBSERVATION futobs;
  double xk[3];
  double **covar;

  ORBFIT_GUARD(guard, result, 3);

  covar = dmatrix(1,6,1,6);

  if (read_abg(abg_file, &p, covar)) orbfit_rethrow();

  /* get observatory code */
  futobs.obscode = obscode;

  futobs.obstime = (jdate - jd0) * DAY;
  futobs.xe = -999.;        /* Force evaluation of earth3d */

  kbo3d_helio(&p, &futobs, xk);

  result[0] = xk[0];
  result[1] = xk[1];
  result[2] = xk[2];

  free_dmatrix(covar,1,6,1,6);
  orbfit_end();
  return result;
}

/* predict.c - Read a file containing a/b/g orbit fit, and spit out
 *  predicted RA & dec plus uncertainties on arbitrary date.
 * 8/12/99 gmb
 */

double *predict(char *abg_file, double jdate, int obscode)
{
  static double result[8];
  jmp_buf guard;
  PBASIS p;
  OBSERVATION	futobs;
  double **covar,**sigxy,a,b,PA,**derivs;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double xx,yy,xy,bovasqrd,det;
  double distance;

  ORBFIT_GUARD(guard, result, 8);

  sigxy = dmatrix(1,2,1,2);
  derivs = dmatrix(1,2,1,2);
  covar = dmatrix(1,6,1,6);
  covecl = dmatrix(1,2,1,2);
  coveq = dmatrix(1,2,1,2);

  if (read_abg(abg_file,&p,covar)) orbfit_rethrow();

  /* get observatory code */
  futobs.obscode=obscode;

  futobs.obstime=(jdate-jd0)*DAY;
  futobs.xe = -999.;		/* Force evaluation of earth3d */

  distance = predict_posn(&p,covar,&futobs,sigxy);

  /* Now transform to RA/DEC, via ecliptic*/
  proj_to_ec(futobs.thetax,futobs.thetay,
	     &lat, &lon,
	     lat0, lon0, derivs);
  /* map the covariance */
  covar_map(sigxy, derivs, covecl, 2, 2);
  
  /* Now to ICRS: */
  ec_to_eq(lat, lon, &ra, &dec, derivs);
  /* map the covariance */
  covar_map(covecl, derivs, coveq, 2, 2);
  
  /* Compute a, b, theta of error ellipse for output */
  xx = coveq[1][1]*cos(dec)*cos(dec);
  xy = coveq[1][2]*cos(dec);
  yy = coveq[2][2];
  PA = 0.5 * atan2(2.*xy,(xx-yy)) * 180./PI;	/*go right to degrees*/
  /* Put PA N through E */
  PA = 90.-PA;
  bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
    / (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
  det = xx*yy-xy*xy;
  b = pow(det*bovasqrd,0.25);
  a = pow(det/bovasqrd,0.25);
  
  ra /= DTOR;
  if (ra<0.) ra+= 360.;
  dec /= DTOR;
  lat /= DTOR;
  lon /= DTOR;
  if (lon<0.) lon+= 360.;

  result[0] = ra;
  result[1] = dec;
  result[2] = a/ARCSEC;
  result[3] = b/ARCSEC;
  result[4] = PA;
  result[5] = distance;
  result[6] = lon;
  result[7] = lat;

  free_dmatrix(sigxy,1,2,1,2);
  free_dmatrix(derivs,1,2,1,2);
  free_dmatrix(covar,1,6,1,6);
  free_dmatrix(covecl,1,2,1,2);
  free_dmatrix(coveq,1,2,1,2);
  orbfit_end();
  return result;
}


double *abg_to_aei(char *abg_file)
{
  static double result[15];
  jmp_buf guard;
  PBASIS p;
  XVBASIS xv;
  ORBIT orbit;
  double d, dd;
  double  **covar_abg, **covar_xyz, **derivs, **covar_aei;

  ORBFIT_GUARD(guard, result, 15);

  covar_abg = dmatrix(1,6,1,6);
  covar_xyz = dmatrix(1,6,1,6);
  covar_aei = dmatrix(1,6,1,6);
  derivs = dmatrix(1,6,1,6);

  if (read_abg(abg_file,&p,covar_abg)) orbfit_rethrow();

  /* Transform the orbit basis and get the deriv. matrix */
  pbasis_to_bary(&p, &xv, derivs);

  /* Map the covariance matrix to new basis */
  covar_map(covar_abg, derivs, covar_xyz,6,6);

  /* Get partial derivative matrix from xyz to aei */
  aei_derivs(&xv, derivs);

  /* Map the covariance matrix to new basis */
  covar_map(covar_xyz, derivs, covar_aei,6,6);

  /* Transform xyz basis to orbital parameters */
  orbitElements(&xv, &orbit);

  d = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
  dd = d*d*sqrt(covar_abg[5][5]);

  result[0] = orbit.a;
  result[1] = orbit.e;
  result[2] = orbit.i;
  result[3] = orbit.lan;
  result[4] = orbit.aop;
  result[5] = orbit.T;
  result[6] = sqrt(covar_aei[1][1]);
  result[7] = sqrt(covar_aei[2][2]);
  result[8] = sqrt(covar_aei[3][3])/DTOR;
  result[9] = sqrt(covar_aei[4][4])/DTOR;
  result[10] = sqrt(covar_aei[5][5])/DTOR;
  result[11] = sqrt(covar_aei[6][6])/DAY;
  result[12] = jd0;
  result[13] = d;
  result[14] = dd;

  free_dmatrix(covar_abg,1,6,1,6);
  free_dmatrix(covar_xyz,1,6,1,6);
  free_dmatrix(covar_aei,1,6,1,6);
  free_dmatrix(derivs,1,6,1,6);
  orbfit_end();
  return result;
}

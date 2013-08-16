#include <Python.h>
/* predict.c - Read a file containing a/b/g orbit fit, and spit out
 *  predicted RA & dec plus uncertainties on arbitrary date.
 * 8/12/99 gmb 
 */
#include "orbfit.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static PyObject *orbfit_fitradec(PyObject *self, PyObject *args)
     //char *mpc_file
{
  const char *mpc_filename;
  const char *abg_filename;
  const char *res_filename;

  FILE *abg_file;
  FILE *mpc_file;
  FILE *res_file;

  OBSERVATION obsarray[MAXOBS];
  int     nobs;


  ORBIT orbit;
  PBASIS p;
  XVBASIS xv;
  
  double **covar;
  double chisq;
  int i; 
  int dof;

  if (!PyArg_ParseTuple(args, "sss", &mpc_filename, &abg_filename, &res_filename ))
    return NULL;

  covar = dmatrix(1,6,1,6);


  if (read_radec(obsarray, mpc_filename, &nobs)) {
    fprintf(stderr, "Error reading input observations\n");
    return NULL ;
  }
  
  /* Call subroutine to do the actual fitting: */
  fit_observations(obsarray, nobs, &p, covar, &chisq, &dof,stdout);

  abg_file = fopen(abg_filename,"w");

  fprintf(abg_file, "# Chi-squared of fit: %.2f DOF: %d\n",chisq,dof);
  fprintf(abg_file, "# Exact a, adot, b, bdot, g, gdot:\n");
  fprintf(abg_file, "%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",p.a,p.adot,p.b,
  	p.bdot, p.g, p.gdot);
  pbasis_to_bary(&p, &xv, NULL);

  orbitElements(&xv, &orbit);
  fprintf(abg_file, "# a=%lf AU,e=%lf,i=%lf deg\n",orbit.a, orbit.e, orbit.i);
  {
    double d, dd;
    d = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
    dd = d*d*sqrt(covar[5][5]);
    fprintf(abg_file, "# Barycentric distance %.3f+-%.3f\n",d,dd);
  }

  /* Print the covariance matrix to stdout */
  fprintf(abg_file, "# Covariance matrix: \n");
  print_matrix(abg_file,covar,6,6);

  /* Print out information on the coordinate system */
  fprintf(abg_file, "#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  fprintf(abg_file, "%12.7f %12.7f %10.7f %10.7f %10.7f  %.6f\n",
	 lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);

  /* Dump residuals to stderr */
  res_file = fopen(res_filename,"w");
  fprintf(res_file,"Best fit orbit gives:\n");
  fprintf(res_file,"obs  time        x      x_resid       y   y_resid\n");
  for (i=0; i<nobs; i++) {
    double x,y;
    kbo2d(&p, &obsarray[i], &x, NULL, &y, NULL);
    fprintf(res_file,"%3d %9.4f %10.3f %7.3f %10.3f %7.3f\n",
	    i, obsarray[i].obstime,
	    obsarray[i].thetax/ARCSEC, (obsarray[i].thetax-x)/ARCSEC,
	    obsarray[i].thetay/ARCSEC, (obsarray[i].thetay-y)/ARCSEC);
  }

  free_dmatrix(covar,1,6,1,6);
  return Py_BuildValue("fffff",orbit.a, orbit.e, orbit.i, orbit.lan, orbit.aop, orbit.T);

} 

  
  

static PyObject *orbfit_predict(PyObject *self, PyObject *args)
     //char *abg_file, float jdate, int obscode)
{
  const char *abg_file;
  float jdate;
  int obscode;

  PBASIS p;
  OBSERVATION	futobs;
  struct date_time dt;
  char	inbuff[256],rastring[20],decstring[20];
  char  outbuff[256];
  char  *f_string;
  double **covar,**sigxy,a,b,PA,**derivs;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double yr,mo,day,hr,mn,ss;
  double xx,yy,xy,bovasqrd,det;
  int i,nfields;
  int iarg=1;

  if (!PyArg_ParseTuple(args,"sfi",&abg_file, &jdate, &obscode))
    return NULL;

  sigxy = dmatrix(1,2,1,2);
  derivs = dmatrix(1,2,1,2);
  covar = dmatrix(1,6,1,6);
  covecl = dmatrix(1,2,1,2);
  coveq = dmatrix(1,2,1,2);

 
  if (read_abg(abg_file,&p,covar) ) { 
    fprintf(stderr, "Error input alpha/beta/gamma file %s\n",abg_file);
    return NULL;
  }

  /* get observatory code */
  futobs.obscode=obscode;

  futobs.obstime=(jdate-jd0)*DAY;
  futobs.xe = -999.;		/* Force evaluation of earth3d */

  predict_posn(&p,covar,&futobs,sigxy);

  /*  This is for the sky plane stuff... */
  /* Compute a, b, theta of error ellipse for output */
  // xx = sigxy[1][1];
  // yy = sigxy[2][2];
  // xy = sigxy[1][2];
  // PA = 0.5 * atan2(2.*xy,(xx-yy)) * 180./PI;	/*go right to degrees*/
  /* Adjust for PA to be N through E, */
  // PA = PA-90;
  // if (PA<-90.) PA += 180.;

  // bovasqrd  = (xx+yy-sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) 
  // (xx+yy+sqrt(pow(xx-yy,2.)+pow(2.*xy,2.))) ;
  // det = xx*yy-xy*xy;
  //b = pow(det*bovasqrd,0.25);
  // a = pow(det/bovasqrd,0.25);

  
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


  return Py_BuildValue("fffff",ra,dec,a/ARCSEC,b/ARCSEC,PA);

}


static PyMethodDef OrbfitMethods[] = {
    {"predict",  orbfit_predict, METH_VARARGS,
     "Given an ABG_FILE, JDATE, OBS_CODE return the objects location."},
    {"fit_radec", orbfit_fitradec, METH_VARARGS,
     "Given an MPC_FILE as input write to ABG_FILE and RES_FILE."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initbk_orbfit(void)
{
  (void) Py_InitModule("bk_orbfit",OrbfitMethods);
}



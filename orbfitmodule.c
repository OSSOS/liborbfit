#include <Python.h>
/* predict.c - Read a file containing a/b/g orbit fit, and spit out
 *  predicted RA & dec plus uncertainties on arbitrary date.
 * 8/12/99 gmb 
 */
#include "orbfit.h"
#include <string.h>
#include <stdlib.h>

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
    return(NULL);
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
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initorbfit(void)
{
  (void) Py_InitModule("orbfit",OrbfitMethods);
}



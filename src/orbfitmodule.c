#include "orbfit.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

double *fitradec(char *mpc_filename, char *abg_filename)
{

  FILE *abg_file ;
  FILE *res_file ;

  OBSERVATION obsarray[MAXOBS];
  int     nobs;

  ORBIT orbit;
  PBASIS p;
  XVBASIS xv;

  double d, dd;
  static double result[2];
  double **covar;
  double chisq;
  int i; 
  int dof;

  covar = dmatrix(1,6,1,6);


  if (read_radec(obsarray, mpc_filename, &nobs)) {
    fprintf(stderr, "Error reading input observations\n");
    return result ;
  }
  
  /* Call subroutine to do the actual fitting: */
  fit_observations(obsarray, nobs, &p, covar, &chisq, &dof, NULL);


  abg_file = fopen(abg_filename,"w");

  /* fprintf(stderr, "# Chi-squared of fit: %.2f DOF: %d\n",chisq,dof); */
  fprintf(abg_file, "# Exact a, adot, b, bdot, g, gdot:\n");
  fprintf(abg_file, "%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",p.a,p.adot,p.b,
  	p.bdot, p.g, p.gdot);
  pbasis_to_bary(&p, &xv, NULL);

  orbitElements(&xv, &orbit);
  /* fprintf(stderr, "# a=%f AU,e=%f,i=%f deg\n",orbit.a, orbit.e, orbit.i); */
  d = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
  dd = d*d*sqrt(covar[5][5]);
  /* fprintf(stderr, "# Barycentric distance %.3f+-%.3f\n",d,dd); */

  /* Print the covariance matrix to the agb_file */
  /* write the covariance to the agbfile */

  fprintf(abg_file, "# Covariance matrix: \n");

  print_matrix(abg_file, covar, 6, 6);



  /* Print out information on the coordinate system */
  fprintf(abg_file, "#     lat0       lon0       xBary     yBary      zBary   JD0\n");
  fprintf(abg_file, "%12.7f %12.7f %10.7f %10.7f %10.7f  %.6f\n",
	 lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);

  fclose(abg_file);

  /* Dump residuals to res_file */
  /*
  res_file = stderr;
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
  */

  free_dmatrix(covar,1,6,1,6);
  result[0] = d;
  result[1] = dd;
  return result;
}


/* predict.c - Read a file containing a/b/g orbit fit, and spit out
 *  predicted RA & dec plus uncertainties on arbitrary date.
 * 8/12/99 gmb
 */

double *predict(char *abg_file, double jdate, int obscode)
{

  PBASIS p;
  OBSERVATION	futobs;
  struct date_time dt;
  char	inbuff[256], rastring[20], decstring[20];
  char  outbuff[256];
  char  *f_string;
  double **covar,**sigxy,a,b,PA,**derivs;
  double lat,lon,**covecl;
  double ra,dec, **coveq;
  double yr,mo,day,hr,mn,ss;
  double xx,yy,xy,bovasqrd,det;
  double distance;
  static double result[6];
  int i,nfields;
  int iarg=1;

  sigxy = dmatrix(1,2,1,2);
  derivs = dmatrix(1,2,1,2);
  covar = dmatrix(1,6,1,6);
  covecl = dmatrix(1,2,1,2);
  coveq = dmatrix(1,2,1,2);

  result[0] = -1.0;
  result[1] = -1.0;
  result[2] = -1.0;
  result[3] = -1.0;
  result[4] = -1.0;
  result[5] = -1.0;
 
  if (read_abg(abg_file,&p,covar) ) { 
    fprintf(stderr, "Error input alpha/beta/gamma file %s\n",abg_file);
    return result;
  }


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


   result[0] = ra;
   result[1] = dec;
   result[2] = a/ARCSEC;
   result[3] = b/ARCSEC;
   result[4] = PA;
   result[5] = distance;

  return result;

}


double *abg_to_aei(char *abg_file)
{
  PBASIS p;
  XVBASIS xv;
  ORBIT orbit;
  static double result[15];
  double d, dd;
  double  **covar_abg, **covar_xyz, **derivs, **covar_aei;

  int	i,j;

  covar_abg = dmatrix(1,6,1,6);
  covar_xyz = dmatrix(1,6,1,6);
  covar_aei = dmatrix(1,6,1,6);
  derivs = dmatrix(1,6,1,6);

  if (read_abg(abg_file,&p,covar_abg)) {
    fprintf(stderr, "Error in input alpha/beta/gamma file\n");
    exit(1);
  }

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

  /* Print out the results, with comments */
  /*
  printf(aei_file,"# Barycentric osculating elements in ICRS at epoch %.1f:\n",jd0);
  printf("#    a            e       i      Node   Arg of Peri   Time of Peri\n");
  printf("%12.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
  	 orbit.a, orbit.e, orbit.i, orbit.lan, orbit.aop, orbit.T);
  fprintf("+-%10.6f  %9.6f  %8.3f %8.3f  %8.3f %11.3f\n",
     sqrt(covar_aei[1][1]),
  	 sqrt(covar_aei[2][2]),
  	 sqrt(covar_aei[3][3])/DTOR,
     sqrt(covar_aei[4][4])/DTOR,
	 sqrt(covar_aei[5][5])/DTOR,
	 sqrt(covar_aei[6][6])/DAY);
  fprintf("# covariance matrix:\n");
  fprint_matrix(stdout,covar_aei,6,6);
  */
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
  return result;
}


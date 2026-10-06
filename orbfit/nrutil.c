#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"

/* Numerical Recipes style vectors and matrices indexed v[nl..nh].  Element 0
 * up to nl-1 is allocated but unused, which avoids forming a pointer before
 * the start of the allocation; only nl >= 0 is supported. */

void nrerror(const char *error_text)
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(int nl, int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	if (nl < 0 || nh < nl) nrerror("bad range in ivector()");
	v=(int *)malloc((size_t) (nh+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v;
}

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	if (nl < 0 || nh < nl) nrerror("bad range in dvector()");
	v=(double *)malloc((size_t) (nh+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	double **m;

	if (nrl < 0 || nrh < nrl || ncl < 0 || nch < ncl)
		nrerror("bad range in dmatrix()");

	/* allocate pointers to rows */
	m=(double **) calloc((size_t) (nrh+1), sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((size_t) (nch+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
	}
	/* return pointer to array of pointers to rows */
	return m;
}

void free_ivector(int *v, int nl, int nh)
/* free an int vector allocated with ivector() */
{
	(void) nl; (void) nh;
	free(v);
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	(void) nl; (void) nh;
	free(v);
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by dmatrix() */
{
	int i;

	(void) ncl; (void) nch;
	if (m == NULL) return;
	for(i=nrh;i>=nrl;i--) free(m[i]);
	free(m);
}

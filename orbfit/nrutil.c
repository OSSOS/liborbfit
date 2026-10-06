#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"
#include "orbfit.h"

/* Numerical Recipes style vectors and matrices indexed v[nl..nh].  Element 0
 * up to nl-1 is allocated but unused, which avoids forming a pointer before
 * the start of the allocation; only nl >= 0 is supported.
 *
 * Every block is kept on a list of live allocations so nr_free_all() can
 * release them after orbfit_fail() unwinds out of a calculation. */

typedef union nrblock {
	struct {
		union nrblock *prev, *next;
	} link;
	max_align_t align;
} nrblock;

static nrblock live = {{&live, &live}};

void *nr_malloc(size_t n)
{
	nrblock *b = malloc(sizeof(nrblock) + n);

	if (!b) nrerror("allocation failure");
	b->link.prev = &live;
	b->link.next = live.link.next;
	live.link.next->link.prev = b;
	live.link.next = b;
	return b + 1;
}

void nr_free(void *p)
{
	nrblock *b;

	if (p == NULL) return;
	b = (nrblock *) p - 1;
	b->link.prev->link.next = b->link.next;
	b->link.next->link.prev = b->link.prev;
	free(b);
}

void nr_free_all(void)
{
	while (live.link.next != &live) nr_free(live.link.next + 1);
}

void nrerror(const char *error_text)
/* Numerical Recipes standard error handler */
{
	orbfit_fail("Numerical Recipes run-time error: %s", error_text);
}

int *ivector(int nl, int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	if (nl < 0 || nh < nl) nrerror("bad range in ivector()");
	return (int *) nr_malloc((size_t) (nh+1)*sizeof(int));
}

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	if (nl < 0 || nh < nl) nrerror("bad range in dvector()");
	return (double *) nr_malloc((size_t) (nh+1)*sizeof(double));
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	double **m;

	if (nrl < 0 || nrh < nrl || ncl < 0 || nch < ncl)
		nrerror("bad range in dmatrix()");

	/* allocate pointers to rows, then the rows */
	m=(double **) nr_malloc((size_t) (nrh+1)*sizeof(double*));
	for(i=0;i<nrl;i++) m[i]=NULL;
	for(i=nrl;i<=nrh;i++)
		m[i]=(double *) nr_malloc((size_t) (nch+1)*sizeof(double));
	return m;
}

void free_ivector(int *v, int nl, int nh)
/* free an int vector allocated with ivector() */
{
	(void) nl; (void) nh;
	nr_free(v);
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	(void) nl; (void) nh;
	nr_free(v);
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by dmatrix() */
{
	int i;

	(void) ncl; (void) nch;
	if (m == NULL) return;
	for(i=nrh;i>=nrl;i--) nr_free(m[i]);
	nr_free(m);
}

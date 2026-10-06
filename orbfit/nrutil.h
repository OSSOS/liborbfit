#ifndef NRUTIL_H
#define NRUTIL_H

double *dvector(int, int);
double **dmatrix(int, int, int, int);
int *ivector(int, int);
void free_dvector(double *, int, int);
void free_ivector(int *, int, int);
void free_dmatrix(double **, int, int, int, int);
void nrerror(const char *);

void lubksb(double **, int, int *, double[]);
void ludcmp(double **a, int n, int *indx, double *d);
void gaussj(double **a, int n, double **b, int m);
void covsrt(double **covar, int ma, int ia[], int mfit);
float ran1(long *idum);
float gasdev(long *idum);

#endif

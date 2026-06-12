#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <R_ext/Lapack.h>

void get_dimnames(int nrow, int ncol,
                  int *irow1, int *irow2, int *irow,
                  int *icol1, int *icol2, int *icol,
                  SEXP dimnames1_, SEXP dimnames2_, SEXP dimnames_);

void get_pos(int nA, int nB, int k, int *i, int *j, int start);

void get_kronecker_index(int nA, int nB, int *iA, int *iB, int ni, int *index, int start);

void sum_set(int n, double *a, double *dx, int *ix, double *b, double *dy, int *iy, double *dz);

void daxty_set(int n, double *a, double *dx, int *ix, double *dy, int *iy, double *dz);

void daxtypz_set(int n, double *a, double *dx, int *ix, double *dy, int *iy, double *dz);

double ddot3_set(int n, double *dx, int *ix, double *dy, int *iy, double *dz);

double dnorm_xty_set(int n, double *dx, int *ix, double *dy, int *iy);

void append_to_order_vector(int k, double *values, int *order);

void get_nu(int N, double *d, int *order, double dmin, double alpha,
            int *nd, double *cumvar, int *nPC, double *totalvar);

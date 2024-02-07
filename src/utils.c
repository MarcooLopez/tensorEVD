#include "tensorEVD.h"

SEXP get_dimnames(int nrow, int ncol,
                  int *irow1, int *irow2, int *irow,
                  int *icol1, int *icol2, int *icol,
                  SEXP dimnames1_,  // optional dimnames from inputs
                  SEXP dimnames2_)
{
  SEXP rownames_ = PROTECT(Rf_allocVector(STRSXP, nrow));
  SEXP colnames_ = PROTECT(Rf_allocVector(STRSXP, ncol));
  SEXP dimnames_ = PROTECT(Rf_allocVector(VECSXP, 2));

  int i, j;
  char *name1 = (char*)malloc(100*sizeof(char));
  char *name2 = (char*)malloc(100*sizeof(char));

  int flag1 = !Rf_isNull(dimnames1_) && (Rf_length(dimnames1_)==2);
  int flag2 = !Rf_isNull(dimnames2_) && (Rf_length(dimnames2_)==2);

  //Rprintf(" Making rownames ...\n");
  for(i=0; i<nrow; i++){
    j = (irow == NULL) ? i : irow[i];
    if(flag1 && (Rf_length(VECTOR_ELT(dimnames1_,0))>0)){
      strcpy(name1, CHAR(STRING_ELT(VECTOR_ELT(dimnames1_,0),irow1[j])));
    }else{
      snprintf(name1, 100, "%d", irow1[j]+1);
    }
    if(flag2 && (Rf_length(VECTOR_ELT(dimnames2_,0))>0)){
      strcpy(name2, CHAR(STRING_ELT(VECTOR_ELT(dimnames2_,0),irow2[j])));
    }else{
      snprintf(name2, 100, "%d", irow2[j]+1);
    }
    SET_STRING_ELT(rownames_,i,mkChar(strcat(strcat(name1, ":"), name2)));
  }

  //Rprintf(" Making colnames ...\n");
  for(i=0; i<ncol; i++){
    j = (icol == NULL) ? i : icol[i];
    if(flag1 && (Rf_length(VECTOR_ELT(dimnames1_,1))>0)){
      strcpy(name1, CHAR(STRING_ELT(VECTOR_ELT(dimnames1_,1),icol1[j])));
    }else{
      snprintf(name1, 100, "%d", icol1[j]+1);
    }
    if(flag2 && (Rf_length(VECTOR_ELT(dimnames2_,1))>0)){
      strcpy(name2, CHAR(STRING_ELT(VECTOR_ELT(dimnames2_,1),icol2[j])));
    }else{
      snprintf(name2, 100, "%d", icol2[j]+1);
    }
    SET_STRING_ELT(colnames_,i,mkChar(strcat(strcat(name1, ":"), name2)));
  }

  SET_VECTOR_ELT(dimnames_, 0, rownames_);
  SET_VECTOR_ELT(dimnames_, 1, colnames_);

  UNPROTECT(3);
  return(dimnames_);
}

//==============================================================
// Obtain the indices for the kronecker product between matrices
// A (nrowA x ncolA) and B (nrowB x ncolB).
// Let n1 and n2 to be either nrow or ncol from A and B,
// the kronecker is obtained by multiplying elements
//     [1,1,...,1,2,2,...,2,...,n1,nA...,nA]  from A
// and elements
//     [1,2,...,nB,1,2,...,nB,...,1,2,...,nB] from B
//
//==============================================================
void get_pos(int nA, int nB, int k, int *i, int *j)
{
  i[0] = floor(k/nB);
  j[0] = k % nB;
}

//==============================================================

void get_kronecker_index(int nA, int nB, int *iA, int *iB, int ni, int *index)
{
  int i, j;

  if(ni == 0){
    int tmp = 0;
    for(i=0; i<nA; i++){
      for(j=0; j<nB; j++){
        iA[tmp] = i;
        iB[tmp] = j;
        tmp++;
      }
    }
  }else{
    for(i=0; i<ni; i++){
      get_pos(nA, nB, index[i], iA+i, iB+i);
    }
  }
}

//====================================================================
// Hadamard product between two vectors:
//       dz[j] <- a * dx[ix[j]] * dy[iy[j]],    j = 1,2,...,n
//
//   [in]  a: (double) A factor to multiply the hadamard by
//   [in]  n: (int) Number of elements in input vector(s) ix and iy
//   [in]  dx: double precision array of dimension <= max(ix)+1
//   [in]  ix: integer array (zero-based) of dimension n
//   [in]  dy: double precision array of dimension <= max(iy)+1
//   [in]  ix: integer array (zero-based) of dimension n
//   [out] dz: double precision array of dimension at least n
//====================================================================
void hadam_set(int n, double *a, double *dx, int *ix, double *dy, int *iy, double *dz)
{
    int m, i;

    /* Clean-up loop so remaining vector length is a multiple of 5.  */
    m = n % 5;
    if(m != 0){
       for(i=0; i<m; i++){
          dz[i] = a[0] * dx[ix[i]] * dy[iy[i]];
       }
       if(n < 5){
          return;
       }
    }
    for(i=m; i<n; i+=5)
    {
       dz[i] = a[0] * dx[ix[i]] * dy[iy[i]];
       dz[i+1] = a[0] * dx[ix[i+1]] * dy[iy[i+1]];
       dz[i+2] = a[0] * dx[ix[i+2]] * dy[iy[i+2]];
       dz[i+3] = a[0] * dx[ix[i+3]] * dy[iy[i+3]];
       dz[i+4] = a[0] * dx[ix[i+4]] * dy[iy[i+4]];
    }
}

//====================================================================
// Euclidean norm of a vector dz that is formed as a Hadamard product
// between two subset vectors:
//       dz[j] <- dx[ix[j]] * dy[iy[j]],    j = 1,2,...,n
// Then the norm is:
//       sqrt(dz[1]^2 + ... + dz[n]^2)
//====================================================================
double dnorm_hadam_set(int n, double *dx, int *ix, double *dy, int *iy)
{
    int m, i;
    double out = 0.0;

    /* Clean-up loop so remaining vector length is a multiple of 5.  */
    m = n % 5;
    if(m != 0){
       for(i=0; i<m; i++){
          out += pow(dx[ix[i]] * dy[iy[i]], 2);
       }
       if(n < 5){
          return(sqrt(out));
       }
    }
    for(i=m; i<n; i+=5)
    {
       out += pow(dx[ix[i]] * dy[iy[i]], 2);
       out += pow(dx[ix[i+1]] * dy[iy[i+1]], 2);
       out += pow(dx[ix[i+2]] * dy[iy[i+2]], 2);
       out += pow(dx[ix[i+3]] * dy[iy[i+3]], 2);
       out += pow(dx[ix[i+4]] * dy[iy[i+4]], 2);
    }

    return(sqrt(out));
}

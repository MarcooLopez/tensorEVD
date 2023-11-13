#include "tensorEVD.h"

//==============================================================
//==============================================================

SEXP R_hadamard(SEXP nrowA_, SEXP ncolA_, SEXP A_,
                SEXP nrowB_, SEXP ncolB_, SEXP B_,
                SEXP nrow_, SEXP ncol_, SEXP out_,
                SEXP irowA_, SEXP icolA_,
                SEXP irowB_, SEXP icolB_,
                SEXP drop_, SEXP kronecker_,
                SEXP makedimnames_)
{
    int nprotect = 5;

    int nrow=INTEGER_VALUE(nrow_);
    int ncol=INTEGER_VALUE(ncol_);
    int nrowA=INTEGER_VALUE(nrowA_);
    int ncolA=INTEGER_VALUE(ncolA_);
    int nrowB=INTEGER_VALUE(nrowB_);
    int ncolB=INTEGER_VALUE(ncolB_);
    int drop=asLogical(drop_);
    int kronecker=asLogical(kronecker_);
    int makedimnames=asLogical(makedimnames_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(B_ = AS_NUMERIC(B_));
    double *B = NUMERIC_POINTER(B_);

    int *irowA, *icolA, *irowB, *icolB;
    if(kronecker)
    { // We take irow/icol indices for the output as those provided in irowA/icolA
      int nirow=Rf_length(irowA_);
      int nicol=Rf_length(icolA_);

      PROTECT(irowA_=AS_INTEGER(irowA_));
      int *irow=INTEGER_POINTER(irowA_);

      PROTECT(icolA_=AS_INTEGER(icolA_));
      int *icol=INTEGER_POINTER(icolA_);

      // Get positions for matrices A and B
      irowA = (int *) R_alloc(nrow, sizeof(int));
      irowB = (int *) R_alloc(nrow, sizeof(int));
      get_kronecker_index(nrowA, nrowB, irowA, irowB, nirow, irow);

      if((nrowA==ncolA) && (nrowB==ncolB) && ((nirow+nicol)==0)){
        // If both are squared matrices, do not repeat the indices
        icolA = irowA;
        icolB = irowB;
      }else{
        icolA=(int *) R_alloc(ncol, sizeof(int));
        icolB=(int *) R_alloc(ncol, sizeof(int));
        get_kronecker_index(ncolA, ncolB, icolA, icolB, nicol, icol);
      }

    }else{
      PROTECT(irowA_=AS_INTEGER(irowA_));
      irowA=INTEGER_POINTER(irowA_);

      PROTECT(irowB_=AS_INTEGER(irowB_));
      irowB=INTEGER_POINTER(irowB_);

      if(Rf_length(icolA_) == 0){
        icolA = irowA;
      }else{
        PROTECT(icolA_=AS_INTEGER(icolA_));
        icolA=INTEGER_POINTER(icolA_);
        nprotect++;
      }

      if(Rf_length(icolB_) == 0){
        icolB = irowB;
      }else{
        PROTECT(icolB_=AS_INTEGER(icolB_));
        icolB=INTEGER_POINTER(icolB_);
        nprotect++;
      }

    }

    // Rprintf(" Making the Hadamard product ...");
    SEXP out2_;
    int ismatrix = 1;
    if((nrow==1) || (ncol==1))
    {
      if(drop){
        out2_ = PROTECT(Rf_allocVector(REALSXP, (long long)nrow*ncol));
        ismatrix = 0;
      }else{
        out2_ = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
      }
    }else{
      out2_ = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
    }
    double *out2 = NUMERIC_POINTER(out2_);

    double a = 1.0;
    size_t j;
    for(j=0; j<ncol; j++){
      hadam_set(nrow, &a, A + (long long)nrowA*icolA[j], irowA, B + (long long)nrowB*icolB[j], irowB, out2 + nrow*j);
    }

    // Rprintf(" Making dimnames ...");
    if(ismatrix && makedimnames){
      setAttrib(out2_, R_DimNamesSymbol, get_dimnames(nrow,ncol,irowA,icolA,irowB,icolB,NULL,NULL));
    }

    UNPROTECT(nprotect);
    return(out2_);
}

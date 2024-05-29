#include "tensorEVD.h"

//==============================================================
// Calculate the entry-wise product
//     A[irowA,icolA]*B[irowB,icolB]
// A is of dimensions: nrowA,ncolA
// B is of dimensions: nrowB,ncolB
//
// For squared matrices (irowA=icolA and irowB=icolB):
//     A[irowA,icolA]*B[irowB,icolB] + C[irowA,icolA]*I[irowB,icolB]
// A and C are of the same dimensions
// B and I are of the same dimensions
//==============================================================
SEXP R_hadamard(SEXP nrowA_, SEXP ncolA_, SEXP A_,
                SEXP nrowB_, SEXP ncolB_, SEXP B_,
                SEXP C_,  // Optional matrix of the same dimension as A
                SEXP irowA_, SEXP icolA_,
                SEXP irowB_, SEXP icolB_,
                SEXP out_, SEXP drop_,
                SEXP makedimnames_, SEXP inplace_)
{
    int nprotect = 4;

    int nrowA = INTEGER_VALUE(nrowA_);
    //int ncolA = INTEGER_VALUE(ncolA_);
    int nrowB = INTEGER_VALUE(nrowB_);
    //int ncolB = INTEGER_VALUE(ncolB_);
    int drop = asLogical(drop_);
    int makedimnames = asLogical(makedimnames_);
    int inplace = INTEGER_VALUE(inplace_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(B_ = AS_NUMERIC(B_));
    double *B = NUMERIC_POINTER(B_);

    int nrow = Rf_length(irowA_);

    PROTECT(irowA_ = AS_INTEGER(irowA_));
    int *irowA = INTEGER_POINTER(irowA_);

    PROTECT(irowB_=AS_INTEGER(irowB_));
    int *irowB = INTEGER_POINTER(irowB_);

    int ncol;
    int *icolA, *icolB;
    if(Rf_length(icolA_) == 0){
      icolA = irowA;
      ncol = nrow;
    }else{
      ncol = Rf_length(icolA_);

      PROTECT(icolA_ = AS_INTEGER(icolA_));
      icolA = INTEGER_POINTER(icolA_);
      nprotect++;
    }

    if(Rf_length(icolB_) == 0){
      icolB = irowB;
    }else{
      PROTECT(icolB_ = AS_INTEGER(icolB_));
      icolB = INTEGER_POINTER(icolB_);
      nprotect++;
    }

    // Rprintf(" Allocating memory for Hadamard product ...\n");
    SEXP out2_;
    int ismatrix = 1;
    double *out2;
    if(inplace == 0){
      // Rprintf(" New memory for a %d x %d matrix\n",nrow,ncol);
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
      out2 = NUMERIC_POINTER(out2_);
      nprotect++;

    }else{
      //out2_ = R_NilValue;
      // Rprintf(" Memory from input %d\n",inplace);
      if(inplace == 1){
        out2 = A;
        out2_ = A_;
      }else{ // inplace==2
        out2 = B;
        out2_ = B_;
      }
    }

    // Rprintf(" Making the Hadamard product ...\n");
    double a = 1.0;
    size_t j;
    for(j=0; j<ncol; j++){
      hadam_set(nrow, &a, A + (long long)nrowA*icolA[j], irowA, B + (long long)nrowB*icolB[j], irowB, out2 + nrow*j);
    }

    if(!Rf_isNull(C_)){
      // Rprintf(" Making shifting ...\n");
      PROTECT(C_ = AS_NUMERIC(C_));
      double *C = NUMERIC_POINTER(C_);
      nprotect++;

      // The output is now a full matrix, thus a common shifting code is applied
      size_t i;
      for(j=0; j<ncol; j++){
        for(i=0; i<nrow; i++){
          if(irowB[i] == icolB[j]){ //if(irowA[i] == icolA[j]){
            //out2[nrow*j + i] += C[nrowB*(long long)icolB[j] + (long long)irowB[i]];
            out2[nrow*j + i] += C[nrowA*(long long)icolA[j] + (long long)irowA[i]];
          }
        }
      }
    }

    if(ismatrix && makedimnames && (inplace==0)){
      // Rprintf(" Making dimnames ...\n");
      SEXP dimnames_ = PROTECT(Rf_allocVector(VECSXP, 2));
      get_dimnames(nrow, ncol, irowA, irowB, NULL, icolA, icolB, NULL,
                   Rf_getAttrib(A_, R_DimNamesSymbol),
                   Rf_getAttrib(B_, R_DimNamesSymbol),
                   dimnames_);
      Rf_setAttrib(out2_, R_DimNamesSymbol, dimnames_);
      //setAttrib(out2_, R_DimNamesSymbol,
      //          get_dimnames(nrow,ncol,irowA,irowB,NULL,icolA,icolB,NULL,
      //                       getAttrib(A_, R_DimNamesSymbol),
      //                       getAttrib(B_, R_DimNamesSymbol)));
      nprotect++;
    }

    UNPROTECT(nprotect);
    return(out2_);
}

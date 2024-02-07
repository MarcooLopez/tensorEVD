#include "tensorEVD.h"

//====================================================================
//  Penalize a numeric covariance matrix formed by the Kronecker
//  product between Sigma and K by adding a penalty formed by the
//  Kronecker between Theta and an identity matrix I.
//  If argument 'byrow' is FALSE (random matrix is stacked by columns)
//            Kronecker(Sigma,K) + Kronecker(Theta,I)
//  Otherwise,
//            Kronecker(K,Sigma) + Kronecker(I,Theta)
//
//====================================================================
SEXP R_kronecker_cov(SEXP nK_, SEXP K_, SEXP nS_, SEXP Sigma_, SEXP Theta_,
                     SEXP byrow_,
                     SEXP nrow_, SEXP ncol_, SEXP out_,
                     SEXP irow_, SEXP icol_, SEXP drop_,
                     SEXP inplace_)
{
    int nprotect = 5;
    int nrow=INTEGER_VALUE(nrow_);
    int ncol=INTEGER_VALUE(ncol_);
    int nK=INTEGER_VALUE(nK_);
    int nS=INTEGER_VALUE(nS_);
    int byrow=asLogical(byrow_);
    int drop=asLogical(drop_);
    int inplace=asLogical(inplace_);
    int nirow=Rf_length(irow_);
    int nicol=Rf_length(icol_);

    PROTECT(K_ = AS_NUMERIC(K_));
    double *K = NUMERIC_POINTER(K_);

    PROTECT(Sigma_ = AS_NUMERIC(Sigma_));
    double *Sigma = NUMERIC_POINTER(Sigma_);

    PROTECT(Theta_=AS_NUMERIC(Theta_));
    double *Theta=NUMERIC_POINTER(Theta_);

    PROTECT(irow_=AS_INTEGER(irow_));
    int *irow=INTEGER_POINTER(irow_);

    PROTECT(icol_=AS_INTEGER(icol_));
    int *icol=INTEGER_POINTER(icol_);

    int *irowS=(int *) R_alloc(nrow, sizeof(int)); // For Sigma
    int *irowK=(int *) R_alloc(nrow, sizeof(int)); // For K
    if(byrow){
      get_kronecker_index(nK, nS, irowK, irowS, nirow, irow);
    }else{
      get_kronecker_index(nS, nK, irowS, irowK, nirow, irow);
    }

    int *icolS, *icolK;
    if((nirow+nicol)==0){ // do not repeat the indices
      icolS = irowS;
      icolK = irowK;
    }else{
      icolS=(int *) R_alloc(ncol, sizeof(int));
      icolK=(int *) R_alloc(ncol, sizeof(int));
      if(byrow){
        get_kronecker_index(nK, nS, icolK, icolS, nicol, icol);
      }else{
        get_kronecker_index(nS, nK, icolS, icolK, nicol, icol);
      }
    }

    SEXP out2_;
    double *out2;
    if(inplace){
      //out2_ = R_NilValue;
      out2_ = K_;
      out2 = K;
    }else{
      if((nrow==1) || (ncol==1))
      {
        if(drop){
          out2_ = PROTECT(Rf_allocVector(REALSXP, (long long)nrow*ncol));
        }else{
          out2_ = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        }
      }else{
        out2_ = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
      }
      out2 = NUMERIC_POINTER(out2_);
      nprotect++;
    }

    //Rprintf(" Making kronecker multiplication by matrix 'a'...\n");
    size_t i, j;
    double a = 1.0;

    for(j=0; j<ncol; j++){
      hadam_set(nrow, &a, K + (long long)nK*icolK[j], irowK, Sigma + (long long)nS*icolS[j], irowS, out2 + nrow*j);
    }

    //Rprintf(" Making shifting ...\n");
    // The output is now a full matrix, thus a common shifting code is applied
    for(j=0; j<ncol; j++){
      for(i=0; i<nrow; i++){
        if(icolK[j] == irowK[i]){
          out2[nrow*j + i] += Theta[nS*(long long)icolS[j] + (long long)irowS[i]];
        }
      }
    }

    UNPROTECT(nprotect);

    return(out2_);
}

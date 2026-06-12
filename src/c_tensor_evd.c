#include "tensorEVD.h"

//==============================================================
//  d1, U1:   Eigenvectors (n1 x nu1) and eigenvalues of K1
//  d2, U2:   Eigenvectors (n2 x nu2) and eigenvalues of K2
//  dmin:     Minimum acepted value for the kronecker eigenvalues
//            (i.e., valuesK1[i]*valuesK2[j])
//  index1:   Which elements from K1 are included in the tensor, zero based
//  index2:   which elements from K2 are included in the tensor, zero based
//  alpha:    Maximum percentage of variance explained
//==============================================================
SEXP R_tensor_evd(SEXP index1_, SEXP index2_,
                  SEXP d1_, SEXP U1_,
                  SEXP d2_, SEXP U2_,
                  SEXP nu_, SEXP alpha_, SEXP dmin_,
                  SEXP norm_, SEXP order_,
                  SEXP makedimnames_, SEXP verbose_)
{
    int nprotect = 10;

    int n = Rf_length(index1_);
    int nu = INTEGER_VALUE(nu_);
    int n1 = Rf_nrows(U1_);
    int nu1 = Rf_ncols(U1_);
    int n2 = Rf_nrows(U2_);
    int nu2 = Rf_ncols(U2_);
    double dmin = NUMERIC_VALUE(dmin_);
    double alpha = NUMERIC_VALUE(alpha_);
    int makedimnames = asLogical(makedimnames_);
    int verbose = asLogical(verbose_);
    long long i, j;

    PROTECT(index1_ = AS_INTEGER(index1_));
    int *index1 = INTEGER_POINTER(index1_);

    PROTECT(index2_ = AS_INTEGER(index2_));
    int *index2 = INTEGER_POINTER(index2_);

    PROTECT(d1_ = AS_NUMERIC(d1_));
    double *d1 = NUMERIC_POINTER(d1_);

    PROTECT(U1_ = AS_NUMERIC(U1_));
    double *U1 = NUMERIC_POINTER(U1_);

    PROTECT(d2_ = AS_NUMERIC(d2_));
    double *d2 = NUMERIC_POINTER(d2_);

    PROTECT(U2_ = AS_NUMERIC(U2_));
    double *U2 = NUMERIC_POINTER(U2_);

    int N = nu1*nu2;
    double *d = (double *) R_alloc(N, sizeof(double));
    int *K1i = (int *) R_alloc(N, sizeof(int)); // Store the Kronecker positions for K1
    int *K2i = (int *) R_alloc(N, sizeof(int)); // Store the Kronecker positions for K2
    double *cumvar = (double *) R_alloc(N, sizeof(double));  // Store accumulated variance

    int flag_norm = Rf_isNull(norm_) ? 0 : 1;
    SEXP norm0_ = PROTECT(Rf_allocVector(REALSXP, N));
    if(flag_norm){
      PROTECT(norm_ = AS_NUMERIC(norm_));
      memcpy(NUMERIC_POINTER(norm0_), NUMERIC_POINTER(norm_), N*sizeof(double));
      nprotect++;
    }
    double *norm = NUMERIC_POINTER(norm0_);

    int flag_order = Rf_isNull(order_) ? 0 : 1;
    SEXP order0_ = PROTECT(Rf_allocVector(INTSXP, N));
    if(flag_order){
      PROTECT(order_ = AS_INTEGER(order_));
      memcpy(INTEGER_POINTER(order0_), INTEGER_POINTER(order_), N*sizeof(int));
      nprotect++;
    }
    int *order = INTEGER_POINTER(order0_);

    // Get the eigenvalues and total variance from the full Kronecker
    // The full Kronecker has nK1*nK2 rows, positions will be saved in K1i and K2i
    // The Kronecker eigenvalues are obtained as the product valuesK1[i]*valuesK2[j]
    // They will be descendantly sorted to return the one with highest variance first
    // until the one that jointly explains certain proportion of total variance
    if(verbose){
      Rprintf(" EVD of K1: n1=%d loadings and nu1=%d eigenvectors\n",n1,nu1);
      Rprintf(" EVD of K2: n2=%d loadings and nu2=%d eigenvectors\n",n2,nu2);
    }

    int k = 0;
    // Loop over the rows of MAP
    if(flag_order)
    { // Order indices are provided
      // The order are shifted to be zero-based
      for(i=0; i<N; i++){
        order[i]--;
      }

      if(flag_norm)
      { // Norms of Hadamard eigenvectors are provided
        if(verbose){
          Rprintf(" Calculating N=%d (nu1 x nu2) eigenvalues ...\n",N);
        }
        for(i=0; i<nu1; i++){
          for(j=0; j<nu2; j++){
            K1i[k] = (int) i;
            K2i[k] = (int) j;

            d[k] = d1[i]*d2[j]*pow(norm[k],2); // normalized Kronecker eigenvalues
            k++;
          }
        }

      }else{ // Obtaining the norms of the Hadamard eigenvectors
        if(verbose){
          Rprintf(" Calculating N=%d (nu1 x nu2) eigenvalues and eigenvectors norms ...\n",N);
        }
        for(i=0; i<nu1; i++){
          for(j=0; j<nu2; j++){
            K1i[k] = (int) i;
            K2i[k] = (int) j;

            norm[k] = dnorm_xty_set(n, U1 + n1*i, index1, U2 + n2*j, index2);
            d[k] = d1[i]*d2[j]*pow(norm[k],2); // normalized Kronecker eigenvalues
            k++;
          }
        }
      }

    }else{  // Use the same procedure but calculating the order using d
      if(flag_norm)
      {
        if(verbose){
          Rprintf(" Calculating N=%d (nu1 x nu2) eigenvalues and orders ...\n",N);
        }
        for(i=0; i<nu1; i++){
          for(j=0; j<nu2; j++){
            K1i[k] = (int) i;
            K2i[k] = (int) j;

            d[k] = d1[i]*d2[j]*pow(norm[k],2); // normalized Kronecker eigenvalues
            append_to_order_vector(k, d, order); // indices of the ordered variances
            k++;
          }
        }

      }else{
        if(verbose){
          Rprintf(" Calculating N=%d (nu1 x nu2) eigenvalues, eigenvectors norms, and orders ...\n",N);
        }
        for(i=0; i<nu1; i++){
          for(j=0; j<nu2; j++){
            K1i[k] = (int) i;
            K2i[k] = (int) j;

            norm[k] = dnorm_xty_set(n, U1 + n1*i, index1, U2 + n2*j, index2);
            d[k] = d1[i]*d2[j]*pow(norm[k],2); // normalized Kronecker eigenvalues

            // Calculate the indices of the ordered variances
            append_to_order_vector(k, d, order);
            k++;
          }
        }
      }
    }

    // Rprintf(" Get the number of PC that accumulated certain variance ...\n");
    double totalvar = 0;
    int nd = N;  // Number of tensor eigenvalues after dropping small ones
    int nu0 = 0; // Number of vectors that accumulates alpha*100% of variance
    get_nu(N, d, order, dmin, alpha, &nd, cumvar, &nu0, &totalvar);

    nu = nu < nu0 ? nu : nu0;

    // Output objects
    SEXP values_ = NULL;
    SEXP vectors_ = NULL;
    if(nu == 0){
      values_ = R_NilValue;
      vectors_ = R_NilValue;

      if(verbose){
        Rprintf(" No tensor eigenvectors were derived ...\n");
      }

    }else{
      if(nd < N){
        if(verbose){
          Rprintf(" Dropped bottom %d of %d eigenvectors with small eigenvalue (< %.5e)\n",N-nd,N,dmin);
        }
      }

      if(verbose){
        Rprintf(" Deriving top nu=%d tensor eigenvectors ...\n",nu);
      }
      // The ith Tensor Eigenvector is formed by the product of the corresponding
      // Eigenvectors of K1 and K2 (positions in K1i and K2i) that were previously
      // sorted by the Tensor variance.
      // Only those entries of the selected Eigenvectors of K1 and K2 that are in
      // the Tensor (provided by indexK1 and indexK2) will enter in the dot product

      // Output objects
      vectors_ = PROTECT(Rf_allocMatrix(REALSXP, n, nu));
      double *vectors = NUMERIC_POINTER(vectors_);

      values_ = PROTECT(Rf_allocVector(REALSXP, nu));
      double *values = NUMERIC_POINTER(values_);
      nprotect+=2;

      double w;
      // Loop over the rows of the (ordered) MAP
      for(i=0; i<nu; i++)
      {
        w = 1/norm[order[i]];
        values[i] = d[order[i]];
        // Derive the ith column Hadamard using selected kronecker eigenvectors and scale it by dividing by w
        daxty_set(n, &w, U1 + (long long)n1*K1i[order[i]], index1, U2 + (long long)n2*K2i[order[i]], index2, vectors + n*i);
      }

      // Set dimnames for vectors
      if(makedimnames){
        SEXP dimnames_ = PROTECT(Rf_allocVector(VECSXP, 2));
        SEXP dimnames1_ = PROTECT(Rf_getAttrib(U1_, R_DimNamesSymbol));
        SEXP dimnames2_ = PROTECT(Rf_getAttrib(U2_, R_DimNamesSymbol));
        get_dimnames(n, nu, index1, index2, NULL, K1i, K2i, order,
                     dimnames1_, dimnames2_, dimnames_);
        Rf_setAttrib(vectors_, R_DimNamesSymbol, dimnames_);
        nprotect+=3;
      }

      if(verbose){
        Rprintf(" Resulting eigenvectors explain %.1f %% of the variance=%f\n",100*cumvar[nu-1],totalvar);
      }
    }

    if(verbose){
      Rprintf(" Done!\n");
    }

    // Return order indices to lie within 1,...,N
    for(i=0; i<N; i++){
      order[i]++;
    }

    SEXP list_ = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(list_, 0, values_);
    SET_VECTOR_ELT(list_, 1, vectors_);
    SET_VECTOR_ELT(list_, 2, norm0_);
    SET_VECTOR_ELT(list_, 3, order0_);
    SET_VECTOR_ELT(list_, 4, ScalarReal(totalvar));

    // Set dimnames for outputs
    SEXP names_ = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(names_, 0, mkChar("values"));
    SET_VECTOR_ELT(names_, 1, mkChar("vectors"));
    SET_VECTOR_ELT(names_, 2, mkChar("norm"));
    SET_VECTOR_ELT(names_, 3, mkChar("order"));
    SET_VECTOR_ELT(names_, 4, mkChar("totalVar"));
    Rf_setAttrib(list_, R_NamesSymbol, names_);

    UNPROTECT(nprotect);

    return(list_);
}

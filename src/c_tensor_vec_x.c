#include "tensorEVD.h"

//==============================================================
//  U1:   Eigenvectors (n1 x nu1) and eigenvalues of K1
//  U2:   Eigenvectors (n2 x nu2) and eigenvalues of K2
//  index1:   Which elements from K1 are included in the tensor, zero based
//  index2:   which elements from K2 are included in the tensor, zero based
//  transpose: If True:  U'(Z1*Z2)' x
//             If False: (Z1*Z2)U x
//==============================================================
SEXP R_tensor_vec_x(SEXP index1_, SEXP index2_,
                    SEXP U1_, SEXP U2_, SEXP x_,
                    SEXP nu_, SEXP order_,
                    SEXP transpose_)
{
    int nprotect = 7;

    int n = Rf_length(index1_);
    int nu = INTEGER_VALUE(nu_);
    int n1 = Rf_nrows(U1_);
    int nu1 = Rf_ncols(U1_);
    int n2 = Rf_nrows(U2_);
    int nu2 = Rf_ncols(U2_);
    int nx = Rf_length(x_);
    int transpose = asLogical(transpose_);
    long long i, j;

    PROTECT(index1_ = AS_INTEGER(index1_));
    int *index1 = INTEGER_POINTER(index1_);

    PROTECT(index2_ = AS_INTEGER(index2_));
    int *index2 = INTEGER_POINTER(index2_);

    PROTECT(x_ = AS_NUMERIC(x_));
    double *x = NUMERIC_POINTER(x_);

    PROTECT(U1_ = AS_NUMERIC(U1_));
    double *U1 = NUMERIC_POINTER(U1_);

    PROTECT(U2_ = AS_NUMERIC(U2_));
    double *U2 = NUMERIC_POINTER(U2_);

    int N = nu1*nu2;
    int *K1i = (int *) R_alloc(N, sizeof(int)); // Store the Kronecker positions for K1
    int *K2i = (int *) R_alloc(N, sizeof(int)); // Store the Kronecker positions for K2
    double *norm = (double *) R_alloc(N, sizeof(double));

    PROTECT(order_ = AS_INTEGER(order_));
    int *order = (int *) R_alloc(N, sizeof(int));
    memcpy(order, INTEGER_POINTER(order_), N*sizeof(int));

    int k = 0;
    // Loop over the rows of MAP
    // Order indices are provided, just calculate the norm

    // The order are shifted to be zero-based
    for(i=0; i<N; i++){
      order[i]--;
    }

    for(i=0; i<nu1; i++)
    {
      for(j=0; j<nu2; j++)
      {
        K1i[k] = (int) i;
        K2i[k] = (int) j;

        // Obtain the norms of the Hadamard eigenvectors
        norm[k] = dnorm_xty_set(n, U1 + n1*i, index1, U2 + n2*j, index2);
        k++;
      }
    }

    if(transpose){
      n = nx < n ? nx : n; // This should not happen
    }else{
      nu = nx < nu ? nx : nu;
    }

    // Output objects
    int n0 = transpose ? nu : n;
    SEXP out_ = PROTECT(Rf_allocVector(REALSXP, n0));
    double *out = NUMERIC_POINTER(out_);

    for(i=0; i<n0; i++){  // Initialize output to zero
      out[i] = 0;
    }

    double w;
    // Loop over the rows of the (ordered) MAP
    if(transpose)
    {
      for(i=0; i<nu; i++){
        w = 1/norm[order[i]];
        out[i] = w * ddot3_set(n, U1 + (long long)n1*K1i[order[i]], index1, U2 + (long long)n2*K2i[order[i]], index2, x);
      }

    }else{
      for(i=0; i<nu; i++){
        w = x[i]/norm[order[i]];
        daxtypz_set(n, &w, U1 + (long long)n1*K1i[order[i]], index1, U2 + (long long)n2*K2i[order[i]], index2, out);
      }
    }

    UNPROTECT(nprotect);
    return(out_);
}

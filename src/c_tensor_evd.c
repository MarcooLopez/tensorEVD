#include "tensorEVD.h"

// ----------------------------------------------------------
// Return sorted (descendant order) indices for a vector of
// numeric values. The ordered position of the element k,
// values[k], is found by comparing with the (ordered) values
// of the previous k-1 values. This is, by finding the last
// index order[j] (j=0,1,...,k) such that
//                values[order[j]] > values[k]
// where 'order' contains the ordered indices for the first
//        k elements of 'values': values[0],...,values[k-1]
// The new vector with sorted indices will contain the k+1
// elements:
//    order[0],...,order[j-1], k, order[j],...,order[k-1]
// ----------------------------------------------------------
void append_to_sorted_vector(int k, double *values, int *order)
{
  int j;

  if(k==0){
    order[0] = k;
  }else{
    j=0;
    while(values[order[j]]>values[k]){
      j++;
      if(j == k){
        break;
      }
    }
    if(j < k){
      memmove(order + j+1, order + j, (k-j)*sizeof(int));
    }
    order[j]=k;
  }
}

//==============================================================
//  VK1:       Eigenvectors of K1
//  VK2:       Eigenvectors of K2
//  valuesK1:  Eigenvalues of K1
//  valuesK2:  Eigenvalues of K2
//  minvalue:  Minimum acepted value for the kronecker eigenvalues
//             (i.e., valuesK1[i]*valuesK2[j])
//  indexK1:   Which elements from K1 are included in the tensor, zero based
//  indexK2:   which elements from K2 are included in the tensor, zero based
//  alpha:     Maximum percentage of variance explained
//==============================================================
SEXP R_tensor_evd(SEXP n_, SEXP n1_, SEXP n2_,
                  SEXP V1_, SEXP V2_,
                  SEXP d1_, SEXP d2_,
                  SEXP minvalue_, SEXP index1_, SEXP index2_,
                  SEXP alpha_,
                  SEXP makedimnames_,
                  SEXP verbose_)
{
    long long i, j;
    double *V2, *V1, *d1, *d2;
    int *index1, *index2;

    int n=INTEGER_VALUE(n_);
    int n1=Rf_length(d1_);
    int n2=Rf_length(d2_);
    double minvalue=NUMERIC_VALUE(minvalue_);
    double alpha=NUMERIC_VALUE(alpha_);
    int makedimnames=asLogical(makedimnames_);
    int verbose=asLogical(verbose_);
    double eps = DBL_EPSILON*100;

    PROTECT(V1_=AS_NUMERIC(V1_));
    V1=NUMERIC_POINTER(V1_);

    PROTECT(V2_=AS_NUMERIC(V2_));
    V2=NUMERIC_POINTER(V2_);

    PROTECT(d1_=AS_NUMERIC(d1_));
    d1=NUMERIC_POINTER(d1_);

    PROTECT(d2_=AS_NUMERIC(d2_));
    d2=NUMERIC_POINTER(d2_);

    PROTECT(index1_=AS_INTEGER(index1_));
    index1=INTEGER_POINTER(index1_);

    PROTECT(index2_=AS_INTEGER(index2_));
    index2=INTEGER_POINTER(index2_);

    int nmap=n1*n2;
    double *d=(double *) R_alloc(nmap, sizeof(double));
    double *w0=(double *) R_alloc(nmap, sizeof(double));
    double *cumvar=(double *) R_alloc(nmap, sizeof(double));
    int *order=(int *) R_alloc(nmap, sizeof(int));
    int *K1i=(int *) R_alloc(nmap, sizeof(int));
    int *K2i=(int *) R_alloc(nmap, sizeof(int));

    // Get the PC's variances and total variance from the full Kronecker
    // The full design Kronecker has nK1*nK2 rows, positions will be saved in K1i and K2i
    // The variance of the tensor PCs are obtained as the product valuesK1[i]*valuesK2[j]
    // They will be descendantly sorted to return the one with highest variance first
    // until the one that jointly explains certain proportion of total variance
    if(verbose){
      Rprintf(" Calculating N=%d (%d x %d) tensor variances ...\n",nmap,n1,n2);
    }
    double totalvar = 0;
    int cont = 0;
    for(i=0; i<n1; i++)  // loop over the rows of MAP
    {
      for(j=0; j<n2; j++)
      {
        K1i[cont] = i;  // Storage the Kronecker positions for K1
        K2i[cont] = j;  // Storage the Kronecker positions for K2

        // Get the norms of the Hadamard eigenvectors
        w0[cont] = dnorm_hadam_set(n, V1 + n1*i, index1, V2 + n2*j, index2);

        // Derive and scale kronecker eigenvalues by multiplying by the squared norm
        d[cont] = d1[i]*d2[j]*pow(w0[cont],2);
    	  totalvar += d[cont];

        // Keep track of the ordered variances
        append_to_sorted_vector(cont, d, order);
        cont++;
      }
    }

    // Rprintf(" Get the number of PC that accumulated certain variance ...\n");
    double cumvar0 = 0; //d[order[0]]/totalvar;
    double mindif = fabs(cumvar0-alpha);
    int nd = nmap;  // N positive tensor eigenvalues
    for(i=0; i<nmap; i++){  // loop over the positive ones to get the minimum
      if(d[order[i]] < minvalue){
        nd = i;
        if(verbose){
          Rprintf(" Dropped bottom %d of %d eigenvectors with eigenvalue smaller than %.5e\n",nmap-nd,nmap,minvalue);
        }
        break;
      }else{
        cumvar[i] = cumvar0 + d[order[i]]/totalvar;
        if(fabs(cumvar[i]-alpha) < mindif){
          mindif = fabs(cumvar[i]-alpha);
        }
        cumvar0 = cumvar[i];
      }
    }
    // Rprintf(" Minimum diff fabs(cumvar[i]-alpha)=%g\n",mindif);

    int nPC = 0;
    for(i=0; i<nd; i++){
      if(fabs(fabs(cumvar[i]-alpha)-mindif) <= eps){
        nPC = i + 1;
        break;
      }
    }

    if(verbose){
      Rprintf(" Top %d of %d eigenvectors explain %.1f %% of the variance=%f\n",nPC,nmap,100*cumvar[nPC-1],totalvar);
    }

    if(verbose){
      Rprintf(" Obtaining tensor eigenvectors ...\n");
    }
    // The ith Tensor Eigenvector is formed by the product of the corresponding
    // Eigenvectors of K1 and K2 (positions in K1i and K2i) that were previously
    // sorted by the Tensor variance.
    // Only those entries of the selected Eigenvectors of K1 and K2 that are in
    // the Tensor (provided by indexK1 and indexK2) will enter in the dot product
    // Output objects
    SEXP vectors_ = PROTECT(Rf_allocMatrix(REALSXP, n, nPC));
    double *vectors=NUMERIC_POINTER(vectors_);

    SEXP values_ = PROTECT(Rf_allocVector(REALSXP, nPC));
    double *values=NUMERIC_POINTER(values_);

    double w;
    for(i=0; i<nPC; i++)  // loop over the rows of the (ordered) MAP
    {
      w = 1/w0[order[i]];
      values[i] = d[order[i]];
      // Derive the ith column Hadamard using selected kronecker eigenvectors and scale it by dividing by w
      hadam_set(n, &w, V1 + (long long)n1*K1i[order[i]], index1, V2 + (long long)n2*K2i[order[i]], index2, vectors + n*i);
    }

    if(verbose){
      Rprintf(" Done!\n");
    }

    // Set dimnames for vectors
    if(makedimnames){
      setAttrib(vectors_, R_DimNamesSymbol,
                get_dimnames(n,nPC,index1,index2,NULL,K1i,K2i,order,
                             getAttrib(V1_, R_DimNamesSymbol),
                             getAttrib(V2_, R_DimNamesSymbol)));
    }

    SEXP list_ = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(list_, 0, values_);
    SET_VECTOR_ELT(list_, 1, vectors_);
    SET_VECTOR_ELT(list_, 2, ScalarReal(totalvar));

    // Set dimnames for outputs
    SEXP names_ = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(names_, 0, mkChar("values"));
    SET_VECTOR_ELT(names_, 1, mkChar("vectors"));
    SET_VECTOR_ELT(names_, 2, mkChar("totalVar"));
    setAttrib(list_, R_NamesSymbol, names_);

    UNPROTECT(10);

    return(list_);
}

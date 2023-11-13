
#====================================================================
# Kronecker penalization:
#
#    Kronecker(Sigma,K) + Kronecker(Theta,I),     if byrow = FALSE
#    Kronecker(K,Sigma) + Kronecker(I,Theta),     otherwise
#====================================================================
# Sigma = NULL; rows <- cols <- NULL; byrow <- FALSE; verbose = TRUE
Kronecker_cov <- function(K, Sigma = NULL, Theta, byrow = FALSE,
                          rows = NULL, cols = NULL, drop = TRUE)
{
  dmK <- dim(K)
  if((sum(dmK)/2)^2 != length(K)){
     stop("'K' must be a squared symmetric matrix")
  }
  nK <- dmK[1]

  dm1 <- dim(Sigma)
  dm2 <- dim(Theta)
  if(length(dm2) == 2L){
    if((sum(dm2)/2)^2 != length(Theta)){
       stop("'Theta' must be a squared symmetric matrix")
    }
  }else{
    if(length(Theta) == 1L){
      if(length(dm1) != 2L){
         stop("'Sigma' must be a squared symmetric matrix")
      }
      Theta <- matrix(as.vector(Theta), nrow=dm1[1], ncol=dm1[2])
    }else{
      Theta <- diag(Theta)
    }
    dm2 <- dim(Theta)
  }
  if(is.null(Sigma)){
    Sigma <- diag(dm2[1])
  }else{
    if((sum(dm1)/2)^2 != length(Sigma)){
       stop("'Sigma' must be a squared symmetric matrix")
    }
  }
  if((nrow(Sigma) != nrow(Theta))|(ncol(Sigma) != ncol(Theta))) {
     stop("'Sigma' and 'Theta' must be of the same dimensions")
  }
  nS <- nrow(Sigma)

  nrows <- nK*nS
  ncols <- nK*nS

  if(!is.null(rows)){
    stopifnot(all(!is.na(rows)))
    if(any(rows<1) | any(rows>nrows)){
      stop("Input 'rows' must be integers between 1 and [nrow(A)*nrow(Sigma)]=",nrows)
    }
    nrows <- length(rows)
  }
  if(!is.null(cols)){
    stopifnot(all(!is.na(cols)))
    if(any(cols<1) | any(cols>ncols)){
      stop("Input 'cols' must be integers between 1 and [ncol(A)*ncol(Sigma)]=",ncols)
    }
    ncols <- length(cols)
  }

  #dyn.load("c_kronecker_cov.so")
  return(.Call('R_kronecker_cov', nK, K, nS, Sigma, Theta, byrow,
                                 nrows, ncols, NULL,
                                 rows-1, cols-1, drop))
  #dyn.unload("c_kronecker_cov.so")

}

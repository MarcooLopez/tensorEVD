
#====================================================================
# Kronecker penalization:
#
#    Kronecker(Sigma,K) + Kronecker(Theta,I),     if byrow = FALSE
#    Kronecker(K,Sigma) + Kronecker(I,Theta),     otherwise
#====================================================================
# Sigma = 1; rows <- cols <- NULL; byrow <- FALSE; verbose = TRUE
Kronecker_cov <- function(K, Sigma = 1, Theta, byrow = FALSE,
                          rows = NULL, cols = NULL, drop = TRUE,
                          inplace = FALSE)
{
  dmK <- dim(K)
  if((sum(dmK)/2)^2 != length(K)){
     stop("'K' must be a squared symmetric matrix")
  }
  nK <- dmK[1]

  scalar_Sigma <- FALSE
  if(length(dim(Sigma)) != 2L){
    if(length(Sigma) == 1L){ # Sigma is an scalar
      scalar_Sigma <- TRUE
      Sigma <- matrix(as.vector(Sigma))
    }else{
      Sigma <- diag(Sigma)
    }
  }
  dm1 <- dim(Sigma)
  if((sum(dm1)/2)^2 != length(Sigma)){
     stop("'Sigma' must be a squared symmetric matrix")
  }

  if(length(dim(Theta)) != 2L){
    if(length(Theta) == 1L){ # Theta is an scalar
      Theta <- matrix(as.vector(Theta), nrow=dm1[1], ncol=dm1[2])
    }else{
      Theta <- diag(Theta)
    }
  }
  dm2 <- dim(Theta)
  if((sum(dm2)/2)^2 != length(Theta)){
     stop("'Theta' must be a squared symmetric matrix")
  }

  if(scalar_Sigma & (dm2[1]>1L)){
    Sigma <- diag(rep(as.vector(Sigma)[1],dm2[1]))
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

  if(inplace){
    if(!(nS==1L) | !(is.null(rows)&is.null(cols))){
      stop("Parameter 'inplace' can be TRUE only when 'Sigma' and 'Theta' are scalars ",
           "and 'rows' and 'cols' are NULL")
    }
  }

  #dyn.load("c_kronecker_cov.so")
  return(.Call('R_kronecker_cov', nK, K, nS, Sigma, Theta, byrow,
                                 nrows, ncols, NULL,
                                 rows-1, cols-1, drop, inplace))
  #dyn.unload("c_kronecker_cov.so")
}

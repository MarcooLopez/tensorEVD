
# EVD1 <- EVD2 <- NULL; d.min = .Machine$double.eps;
# make.dimnames = FALSE; normalize=verbose = TRUE
tensorEVD <- function(K1, K2, ID1, ID2, EVD1 = NULL,
                      EVD2 = NULL, nu = NULL, alpha = 1,
                      d.min = .Machine$double.eps*10,
                      norm = NULL, order = NULL,
                      make.dimnames = FALSE, verbose = FALSE)
{
    isEigen1 <- isEigen2 <- FALSE

    stopifnot(length(ID1) == length(ID2))

    if(alpha<0 | alpha>1){
      stop("Parameter 'alpha' must be a number between 0 and 1")
    }

    flag1 <- ifelse(is.list(EVD1),all(c("values","vectors")%in%names(EVD1))&(length(EVD1)==2),FALSE)
    flag2 <- ifelse(is.list(EVD2),all(c("values","vectors")%in%names(EVD2))&(length(EVD2)==2),FALSE)

    # For K1
    if(flag1){
      if(!missing(K1)){
        message(" Both 'K1' and 'EVD1' are provided. 'EVD1' will be used and 'K1' will be ignored")
      }
      dm1 <- dim(EVD1$vectors)
      if(length(dm1) != 2L){
         stop("'EVD1$vectors' must be a matrix with eigenvectors in columns")
      }
      if(dm1[2] != length(EVD1$values)){
        stop("'ncol(EVD1$vectors)' must be the same as 'length(EVD1$values)'")
      }
      index1 <- match_ID(EVD1$vectors, ID1, check=FALSE)
      isEigen1 <- TRUE

    }else{
      index1 <- match_ID(K1, ID1, check=TRUE)

      EVD1 <- eigen(K1, symmetric=TRUE)
      if(has_names(K1)){
        rownames(EVD1$vectors) <- rownames(K1)
      }
      dm1 <- dim(EVD1$vectors)
    }

    # For K2
    if(flag2){
      if(!missing(K2)){
        message(" Both 'K2' and 'EVD2' are provided. 'EVD2' will be used and 'K2' will be ignored")
      }
      dm2 <- dim(EVD2$vectors)
      if(length(dm2) != 2L){
         stop("'EVD2$vectors' must be a matrix with eigenvectors in columns")
      }
      if(dm2[2] != length(EVD2$values)){
        stop("'ncol(EVD2$vectors)' must be the same as 'length(EVD2$values)'")
      }
      index2 <- match_ID(EVD2$vectors, ID2, check=FALSE)
      isEigen2 <- TRUE

    }else{
      index2 <- match_ID(K2, ID2, check=TRUE)

      EVD2 <- eigen(K2, symmetric=TRUE)
      if(has_names(K2)){
        rownames(EVD2$vectors) <- rownames(K2)
      }
      dm2 <- dim(EVD2$vectors)
    }

    if(is.null(index1)){
      stop("'ID1' could not be matched to ",ifelse(isEigen1,"rows of 'EVD1$vectors'","rows/columns of 'K1'"))
    }
    if(is.null(index2)){
      stop("'ID2' could not be matched to ",ifelse(isEigen1,"rows of 'EVD2$vectors'","rows/columns of 'K2'"))
    }

    N <- dm1[2]*dm2[2]  # Size of the Kronecker

    if(!is.null(norm)){
      if(length(norm) != N){
        stop("'norm' should be a vector of length equal to ncol(EVD1$vectors)*ncol(EVD2$vectors) = ",N)
      }
    }

    if(!is.null(order)){
      flag <- (min(order) == 1) & (max(order) == N)
      if((length(order) != N) | !flag){
        stop("'order' should be a vector of length equal to ncol(EVD1$vectors)*ncol(EVD2$vectors) = ",N,",\n",
             " with permutation indices 1,2,...,",N)
      }
    }

    nu <- ifelse(is.null(nu), N, nu)
    if(nu<0 | nu>N){
      stop("'nu' must be an integer between 0 and ",N)
    }

    #dyn.load("c_tensor_evd.so")
    return(.Call('R_tensor_evd', index1, index2,
                                 EVD1$values, EVD1$vectors,
                                 EVD2$values, EVD2$vectors,
                                 nu, as.numeric(alpha), as.numeric(d.min),
                                 norm, order, make.dimnames, verbose))
    #dyn.unload("c_tensor_evd.so")
}

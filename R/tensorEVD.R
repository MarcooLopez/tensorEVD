
# V1 <- d1 <- V2 <- d2 <- NULL; d.min = .Machine$double.eps; verbose= TRUE
tensorEVD <- function(K1, K2, ID1, ID2, alpha = 1.0,
                      V1 = NULL, d1 = NULL, V2 = NULL, d2 = NULL,
                      d.min = .Machine$double.eps,
                      make.dimnames = FALSE,
                      verbose = FALSE){

    isEigen1 <- isEigen2 <- FALSE
    names1 <- names2 <- NULL
    # For K1
    if(missing(K1)){
      if(is.null(V1) | is.null(d1)){
        stop("'V1' and 'd1' must be provided when 'K1' is missing")
      }
      dm1 <- dim(V1)
      if((sum(dm1)/2)^2 != length(V1)){
         stop("'V1' must be a squared matrix with eigenvectors in columns")
      }
      stopifnot(dm1[2] == length(d1))
      names1 <- rownames(V1)
      isEigen1 <- TRUE

    }else{
      dm1 <- dim(K1)
      if((sum(dm1)/2)^2 != length(K1)){
         stop("'K1' must be a squared symmetric matrix")
      }

      names1 <- rownames(K1)
      if(has_names(K1)){
        stopifnot(all(rownames(K1) == colnames(K1)))
      }
    }

    # For K2
    if(missing(K2)){
      if(is.null(V2) | is.null(d2)){
        stop("'V2' and 'd2' must be provided when 'K2' is missing")
      }
      dm2 <- dim(V2)
      if((sum(dm2)/2)^2 != length(V2)){
         stop("'V2' must be a squared matrix with eigenvectors in columns")
      }
      stopifnot(dm2[2] == length(d2))
      names2 <- rownames(V2)
      isEigen2 <- TRUE

    }else{
      dm2 <- dim(K2)
      if((sum(dm2)/2)^2 != length(K2)){
         stop("'K2' must be a squared symmetric matrix")
      }

      names2 <- rownames(K2)
      if(has_names(K2)){
        stopifnot(all(rownames(K2) == colnames(K2)))
      }
    }

    stopifnot(length(ID1) == length(ID2))
    n <- length(ID1)

    index1 <- get_index(dm1[1], names1, ID1)
    index2 <- get_index(dm2[1], names2, ID2)

    if(is.null(index1)){
      stop("'ID1' could not be matched to 'K1'")
    }
    if(is.null(index2)){
      stop("'ID2' could not be matched to 'K2'")
    }

    if(isEigen1){
      EVD1 <- list(vectors=V1, values=d1)
    }else{
      EVD1 <- eigen(K1, symmetric=TRUE)
    }

    if(isEigen2){
      EVD2 <- list(vectors=V2, values=d2)
    }else{
      EVD2 <- eigen(K2, symmetric=TRUE)
    }

    #dyn.load("c_tensor_evd.so")
    return(.Call('R_tensor_evd', n, dm1[1], dm2[1],
                                 EVD1$vectors, EVD2$vectors, EVD1$values, EVD2$values,
                                 as.numeric(d.min), index1-1L, index2-1L,
                                 as.numeric(alpha), make.dimnames, verbose))
    #dyn.unload("c_tensor_evd.so")

}

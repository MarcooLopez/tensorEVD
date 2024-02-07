
# EVD1 <- EVD2 <- NULL; d.min = .Machine$double.eps; verbose= TRUE
tensorEVD <- function(K1, K2, ID1, ID2, alpha = 1.0,
                      EVD1 = NULL, EVD2 = NULL,
                      d.min = .Machine$double.eps,
                      make.dimnames = FALSE,
                      verbose = FALSE)
{
    isEigen1 <- isEigen2 <- FALSE
    names1 <- names2 <- NULL
    # For K1
    if(missing(K1)){
      flag <- ifelse(is.list(EVD1),all(c("values","vectors")%in%names(EVD1))&(length(EVD1)==2),FALSE)
      if(!flag){
        stop("A list type object 'EVD1' must be provided when 'K1' is missing.",
             "\nThis should contain 'values' and 'vectors' as per the 'eigen' function")
      }
      dm1 <- dim(EVD1$vectors)
      if((sum(dm1)/2)^2 != length(EVD1$vectors)){
         stop("'EVD1$vectors' must be a squared matrix with eigenvectors in columns")
      }
      stopifnot(dm1[2] == length(EVD1$values))
      names1 <- rownames(EVD1$vectors)
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
      flag <- ifelse(is.list(EVD2),all(c("values","vectors")%in%names(EVD2))&(length(EVD2)==2),FALSE)
      if(!flag){
        stop("A list type object 'EVD2' must be provided when 'K2' is missing.",
             "\nThis should contain 'values' and 'vectors' as per the 'eigen' function")
      }
      dm2 <- dim(EVD2$vectors)
      if((sum(dm2)/2)^2 != length(EVD2$vectors)){
         stop("'EVD2$vectors' must be a squared matrix with eigenvectors in columns")
      }
      stopifnot(dm2[2] == length(EVD2$values))
      names2 <- rownames(EVD2$vectors)
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

    if(!isEigen1){
      EVD1 <- eigen(K1, symmetric=TRUE)
      if(has_names(K1)){
        rownames(EVD1$vectors) <- rownames(K1)
      }
    }
    if(!isEigen2){
      EVD2 <- eigen(K2, symmetric=TRUE)
      if(has_names(K2)){
        rownames(EVD2$vectors) <- rownames(K2)
      }
    }

    #dyn.load("c_tensor_evd.so")
    return(.Call('R_tensor_evd', n, dm1[1], dm2[1],
                                 EVD1$vectors, EVD2$vectors, EVD1$values, EVD2$values,
                                 as.numeric(d.min), index1-1L, index2-1L,
                                 as.numeric(alpha), make.dimnames, verbose))
    #dyn.unload("c_tensor_evd.so")

}

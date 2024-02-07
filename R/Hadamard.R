
#====================================================================
# Hadamard product between matrices A and B
#====================================================================
# colsA = NULL; colsB = NULL; drop <- TRUE; make.dimnames <- inplace <- FALSE
Hadamard <- function(A, B, rowsA, rowsB,
                     colsA = NULL, colsB = NULL,
                     make.dimnames = FALSE, drop = TRUE,
                     inplace = FALSE)
{

  if((length(dim(A)) != 2L)){
    A <- as.matrix(A, ncol=1L)
  }
  if((length(dim(B)) != 2L)){
    B <- as.matrix(B, ncol=1L)
  }

  dmA <- dim(A)
  dmB <- dim(B)

  # Match rows IDs
  fixedA <- fixedB <- c(FALSE,FALSE)
  if(missing(rowsA)){
    indexrowA <- seq(0,dmA[1]-1)
    fixedA[1] <- TRUE
  }else{
    indexrowA <- get_index(dmA[1], rownames(A), rowsA)
    if(is.null(indexrowA)){
      stop("'rowsA' could not be matched to rows of 'A'")
    }
    indexrowA <- indexrowA-1  # to pass zero-based indices
  }

  if(missing(rowsB)){
    indexrowB <- seq(0,dmB[1]-1)
    fixedB[1] <- TRUE
  }else{
    indexrowB <- get_index(dmB[1], rownames(B), rowsB)
    if(is.null(indexrowB)){
      stop("'rowsB' could not be matched to rows of 'B'")
    }
    indexrowB <- indexrowB-1  # to pass zero-based indices
  }

  # Checkpoint for rows IDs
  if(length(indexrowA) != length(indexrowB)){
    tmp <- c(ifelse(fixedA[1],"'rowsA'",NA),ifelse(fixedB[1],"'rowsB'",NA))
    tmp <- paste0(".\nOtherwise, provide parameter(s) ",paste(tmp[!is.na(tmp)],collapse=" and "))
    tmp <- ifelse(any(c(fixedA[1],fixedB[1])),tmp,"")
    stop(ifelse(fixedA[1],"'nrow(A)'","'length(rowsA)'")," should be the same as ",
         ifelse(fixedB[1],"'nrow(B)'","'length(rowsB)'"),tmp)
  }
  nrows <- length(indexrowA)

  # Match columns IDs
  indexcolA <- indexcolB <- NULL
  if(is.null(colsA)){
    tmp <- (dmA[1]==dmA[2]) & has_names(A)
    flag <- ifelse(tmp, all(rownames(A)==colnames(A)), FALSE)
    if(flag){
      indexcolA <- indexrowA
      fixedA[2] <- fixedA[1]
    }else{
      indexcolA <- seq(0,dmA[2]-1)
      fixedA[2] <- TRUE
    }
  }else{
    indexcolA <- get_index(dmA[2], colnames(A), colsA)
    if(is.null(indexcolA)){
      stop("'colsA' could not be matched to columns of 'A'")
    }
    indexcolA <- indexcolA-1  # to pass zero-based indices
  }

  if(is.null(colsB)){
    tmp <- (dmB[1]==dmB[2]) & has_names(B)
    flag <- ifelse(tmp, all(rownames(B)==colnames(B)), FALSE)
    if(flag){
      indexcolB <- indexrowB
      fixedB[2] <- fixedB[1]
    }else{
      indexcolB <- seq(0,dmB[2]-1)
      fixedB[2] <- TRUE
    }
  }else{
    indexcolB <- get_index(dmB[2], colnames(B), colsB)
    if(is.null(indexcolB)){
      stop("'colsB' could not be matched to columns of 'B'")
    }
    indexcolB <- indexcolB-1  # to pass zero-based indices
  }

  if(length(indexcolA) != length(indexcolB)){
    tmp <- c(ifelse(is.null(colsA)&fixedA[2],"'colsA'",NA),ifelse(is.null(colsB)&fixedB[2],"'colsB'",NA))
    tmp <- paste0(".\nOtherwise, provide parameter(s) ",paste(tmp[!is.na(tmp)],collapse=" and "))
    tmp <- ifelse((is.null(colsA)&fixedA[2])|(is.null(colsB)&fixedB[2]),tmp,"")
    tt1 <- ifelse(is.null(colsA),ifelse(fixedA[2],"ncol(A)","length(rowsA)"),"length(colsA)")
    tt2 <- ifelse(is.null(colsB),ifelse(fixedB[2],"ncol(B)","length(rowsB)"),"length(colsB)")
    stop("'",tt1,"' should be the same as '",tt2,"'",tmp)
  }
  ncols <- length(indexcolA)

  if(inplace){
    inplace2 <- ifelse(all(fixedA),1,ifelse(all(fixedB),2,0))
    if(inplace2 == 0){
      stop("'inplace' calculation can be only applied when either 'A' or 'B' are not resized as per ",
           "the 'rows' and 'cols' parameters")
    }
  }else{
    inplace2 <- 0
  }

  #dyn.load("c_hadamard.so")
  return(.Call('R_hadamard', dmA[1], dmA[2], A, dmB[1], dmB[2], B,
                             nrows, ncols, NULL,
                             indexrowA, indexcolA,
                             indexrowB, indexcolB,
                             drop, FALSE, make.dimnames, inplace2))
  #dyn.unload("c_hadamard.so")

}

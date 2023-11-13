
#====================================================================
# Hadamard product between matrices A and B
#====================================================================
# colsA = NULL; colsB = NULL; drop <- unpack <- TRUE; verbose=TRUE
Hadamard <- function(A, B, rowsA, rowsB,
                     colsA = NULL, colsB = NULL,
                     make.dimnames = FALSE, drop = TRUE)
{

  if((length(dim(A)) != 2L)){
    A <- matrix(A, ncol=1L)
  }
  if((length(dim(B)) != 2L)){
    B <- matrix(B, ncol=1L)
  }

  dmA <- dim(A)
  dmB <- dim(B)

  if(is.null(colsA)){
    if(dmA[1] != dmA[2]){
      stop("'colsA' must be provided when 'A' is not a squared matrix")
    }
    if(has_names(A)){
      stopifnot(all(rownames(A) == colnames(A)))
    }
  }

  if(is.null(colsB)){
    if(dmB[1] != dmB[2]){
      stop("'colsB' must be provided when 'B' is not a squared matrix")
    }
    if(has_names(B)){
      stopifnot(all(rownames(B) == colnames(B)))
    }
  }

  # Match rows IDs
  stopifnot(length(rowsA) == length(rowsB))
  nrows <- length(rowsA)

  indexrowA <- get_index(dmA[1], rownames(A), rowsA)
  indexrowB <- get_index(dmB[1], rownames(B), rowsB)

  if(is.null(indexrowA)){
    stop("'rowsA' could not be matched to rows of 'A'")
  }
  if(is.null(indexrowB)){
    stop("'rowsB' could not be matched to rows of 'B'")
  }

  # Match columns IDs
  if(!is.null(colsA) & !is.null(colsB)){
    stopifnot(length(colsA) == length(colsB))
    ncols <- length(colsA)
  }else{
    if(!is.null(colsA)){
      stopifnot(length(colsA) == nrows)
    }
    if(!is.null(colsB)){
      stopifnot(length(colsB) == nrows)
    }
    ncols <- nrows
  }

  indexcolA <- indexcolB <- NULL
  if(!is.null(colsA)){
    indexcolA <- get_index(dmA[2], colnames(A), colsA)
    if(is.null(indexcolA)){
      stop("'colsA' could not be matched to columns of 'A'")
    }
  }
  if(!is.null(colsB)){
    indexcolB <- get_index(dmB[2], colnames(B), colsB)
    if(is.null(indexcolB)){
      stop("'colsB' could not be matched to columns of 'B'")
    }
  }

  #dyn.load("c_hadamard.so")
  return(.Call('R_hadamard', dmA[1], dmA[2], A, dmB[1], dmB[2], B,
                             nrows, ncols, NULL,
                             indexrowA-1, indexcolA-1,
                             indexrowB-1, indexcolB-1,
                             drop, FALSE, make.dimnames))
  #dyn.unload("c_hadamard.so")

}

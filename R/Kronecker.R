
#====================================================================
# Kronecker product between matrices A and B
#====================================================================
# rows <- cols <- NULL; drop = TRUE; unpack = TRUE; verbose=TRUE
Kronecker <- function(A, B, rows = NULL, cols = NULL,
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
  nrows <- dmA[1]*dmB[1]
  ncols <- dmA[2]*dmB[2]

  if(!is.null(rows)){
    stopifnot(all(!is.na(rows)))
    if(any(rows<1) | any(rows>nrows)){
      stop("Input 'rows' must be integers between 1 and nrow(A)*nrow(B)=",nrows)
    }
    nrows <- length(rows)
  }
  if(!is.null(cols)){
    stopifnot(all(!is.na(cols)))
    if(any(cols<1) | any(cols>ncols)){
      stop("Input 'cols' must be integers between 1 and ncol(A)*ncol(B)=",ncols)
    }
    ncols <- length(cols)
  }

  #dyn.load("c_hadamard.so")
  return(.Call('R_hadamard', dmA[1], dmA[2], A, dmB[1], dmB[2], B,
                             nrows, ncols, NULL,
                             rows-1, cols-1, NULL, NULL,
                             drop, TRUE, make.dimnames))
  #dyn.unload("c_hadamard.so")

}

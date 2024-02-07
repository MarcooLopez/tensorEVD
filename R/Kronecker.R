
#====================================================================
# Kronecker product between matrices A and B
#====================================================================
# rows <- cols <- NULL; drop = TRUE; inplace=TRUE
Kronecker <- function(A, B, rows = NULL, cols = NULL,
                      make.dimnames = FALSE, drop = TRUE,
                      inplace = FALSE)
{

  if(length(dim(A)) != 2L){
    A <- as.matrix(A, ncol=1L)
  }
  if(length(dim(B)) != 2L){
    B <- as.matrix(B, ncol=1L)
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
    rows <- rows-1
  }
  if(!is.null(cols)){
    stopifnot(all(!is.na(cols)))
    if(any(cols<1) | any(cols>ncols)){
      stop("Input 'cols' must be integers between 1 and ncol(A)*ncol(B)=",ncols)
    }
    ncols <- length(cols)
    cols <- cols-1
  }

  if(inplace){
    inplace2 <- ifelse((dmB[1]*dmB[2])==1,1,ifelse((dmA[1]*dmA[2])==1,2,0))
    if(!(inplace2>0) | !(is.null(rows)&is.null(cols))){
      stop("'inplace' calculation can be only applied when either 'A' or 'B' are not resized:",
           "\n one of them is a scalar, and 'rows' and 'cols' are NULL")
    }
  }else{
    inplace2 <- 0
  }

  #dyn.load("c_hadamard.so")
  return(.Call('R_hadamard', dmA[1], dmA[2], A, dmB[1], dmB[2], B,
                             nrows, ncols, NULL,
                             rows, cols, NULL, NULL,
                             drop, TRUE, make.dimnames, inplace2))
  #dyn.unload("c_hadamard.so")

}

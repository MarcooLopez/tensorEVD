
# U1 <- d1 <- U2 <- d2 <- NULL; verbose= TRUE
tensorVec_x <- function(U1, U2, ID1, ID2, x, order,
                        nu = NULL, transpose = FALSE)
{

    stopifnot(length(ID1) == length(ID2))

    # For K1
    dm1 <- dim(U1)
    if(length(dm1) != 2L){
       stop("'U1' must be a matrix with eigenvectors in columns")
    }
    index1 <- match_ID(U1, ID1, check=FALSE)

    # For K2
    dm2 <- dim(U2)
    if(length(dm2) != 2L){
      stop("'U2' must be a matrix with eigenvectors in columns")
    }
    index2 <- match_ID(U2, ID2, check=FALSE)

    if(is.null(index1)){
      stop("'ID1' could not be matched to rows of 'U1'")
    }
    if(is.null(index2)){
      stop("'ID2' could not be matched to rows of 'U2'")
    }

    N <- dm1[2]*dm2[2]  # Size of the Kronecker
    if(!is.null(order)){
      flag <- (min(order) == 1) & (max(order) == N)
      if((length(order) != N) | !flag){
        stop("'order' should be a vector of length equal to ncol(U1)*ncol(U2) = ",N,",\n",
             " with permutation indices 1,2,...,",N)
      }
    }

    if(transpose){
      if(length(x) != length(ID1)){
        stop("'x' should be a numeric vector of length equal to length(ID1) = ",length(ID1))
      }
    }else{
      if(length(x) > N){
        stop("'x' should be a numeric vector of length at most ncol(U1)*ncol(U2) = ",N)
      }
    }

    nu <- ifelse(is.null(nu), N, nu)
    if(nu<0 | nu>N){
      stop("'nu' must be an integer between 0 and ",N)
    }

    #dyn.load("c_tensor_vec_x.so")
    return(.Call('R_tensor_vec_x', index1, index2,
                                   U1, U2, x, nu,
                                   order, transpose))
    #dyn.unload("c_tensor_vec_x.so")
}

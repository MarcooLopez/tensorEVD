
#====================================================================
#====================================================================

has_names <- function(A){
  dm <- dim(A)
  if(length(dm) == 2L){
    out <- length(unlist(dimnames(A))) == sum(dm)
  }else{
    out <- FALSE
  }

  out
}

#====================================================================
#====================================================================

get_index <- function(n, Names, ID){
  index <- NULL
  stopifnot(all(!is.na(ID)))

  if(is.null(Names)){
    if(is.numeric(ID) | is.integer(ID)){
      ID <- as.integer(ID)
      rg <- range(ID)
      if((1L <= rg[1]) & (rg[2] <= n)){
        index <- ID
      }
    }
  }else{
    if(all(as.character(ID) %in% Names)){
      index <- match(as.character(ID), Names)
    }else{
      if(is.numeric(ID) | is.integer(ID)){
        ID <- as.integer(ID)
        rg <- range(ID)
        if((1L <= rg[1]) & (rg[2] <= n)){
          index <- ID
        }
      }
    }
  }
  return(index)
}

#====================================================================
#====================================================================

capitalize <- function(string){
  substr(string,1,1) <- toupper(substr(string,1,1))
  string
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
  |=======================================================================|
  | Loaded '",pkgname,"' R-package. Version ", utils::packageVersion(pkgname),"
  | Authors: Lopez-Cruz M, Perez-Rodriguez P, & de los Campos G
  |=======================================================================|
  ")

  tmp <- utils::old.packages(repos="https://cloud.r-project.org")
  if(pkgname %in% rownames(tmp)){
    packageStartupMessage(" Note: New version ",tmp[pkgname,"ReposVer"],
            " of this package is available on CRAN")
  }
}


\name{Multivariate variance matrix}
\alias{Kronecker_cov}
\alias{Hadamard_cov}
\title{Multivariate variance matrix penalization}
\usage{
Kronecker_cov(Sigma = 1, K, Theta, swap = FALSE,
              rows = NULL, cols = NULL, 
              drop = TRUE, inplace = FALSE)
              
Hadamard_cov(Sigma = 1, K, Theta, IDS, IDK, 
             drop = TRUE, inplace = FALSE) 
              
}
\arguments{
\item{Sigma}{(numeric) A variance matrix among features. If is scalar, a scaled identity matrix with the same dimension as \code{Theta} is used}

\item{K}{(numeric) Variance matrix among subjects}

\item{Theta}{(numeric) A diagonal-shifting parameter, value to be added to the diagonals of the resulting (co)variance matrix. It should be a (symmetric) matrix with the same dimension as \code{Sigma}}

\item{rows}{(integer) Index which rows of the (Kronecker product) (co)variance matrix are to be returned. Default \code{rows=NULL} will return all the rows}

\item{cols}{(integer) Index which columns of the (Kronecker product) (co)variance are to be returned. Default \code{cols=NULL} return all the columns}

\item{IDS}{(integer/character) Vector with either indices or row names mapping from rows/columns of \code{Sigma} and \code{Theta} into the resulting (Hadamard product) (co)variance matrix}

\item{IDK}{(integer/character) Vector with either indices or row names mapping from rows/columns of \code{K} into the resulting (Hadamard product) (co)variance matrix}

\item{swap}{(logical) Either \code{TRUE} or \code{FALSE} (default) to whether swap the order of the matrices in the resulting (Kronecker product) (co)variance matrix}

\item{drop}{(logical) Either \code{TRUE} or \code{FALSE} to whether return a uni-dimensional vector when output is a matrix with either 1 row or 1 column as per the \code{rows} and \code{cols} arguments}

\item{inplace}{(logical) Either \code{TRUE} or \code{FALSE} to whether operate directly on matrix \code{K} when \code{Sigma} and \code{Theta} are scalars. This is possible only when \code{rows=NULL} and \code{cols=NULL}. When \code{TRUE} the output will be overwritten on the same address occupied by \code{K}. Default \code{inplace=FALSE}}

}
\value{
Returns the penalized (co)variance matrix formed either as a Kronecker or Hadamard product. For the Kronecker product case, it can be a sub-matrix of the Kronecker product as per the \code{rows} and \code{cols} arguments.
}
\description{
Ridge penalization of a multi-variate (co)variance matrix taking the form of either a Kronecker or Hadamard product
}
\details{
Assume that a multi-variate random matrix \ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}} with 
\ifelse{html}{\out{<i>n</i>}}{\eqn{n}{n}} subjects in rows and \ifelse{html}{\out{<i>p</i>}}{\eqn{p}{p}} features in columns follows a matrix Gaussian distribution with certain matrix of means 
\ifelse{html}{\out{<b>M</b>}}{\eqn{\textbf{M}}{M}} and variance matrix
\ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} of dimension 
\ifelse{html}{\out{<i>n</i> &times; <i>n</i>}}{\eqn{n\times n}{n x n}} between subjects, and 
\ifelse{html}{\out{<b>&Sigma;</b>}}{\eqn{\boldsymbol{\Sigma}}{Sigma}} of dimension 
\ifelse{html}{\out{<i>p</i> &times; <i>p</i>}}{\eqn{p\times p}{p x p}} between features.

\ifelse{html}{\out{<b>Kronecker product form.</b>}}{\eqn{\textbf{Kronecker product form.}}{Kronecker product form.}}

The random variable 
\ifelse{html}{\out{<b>x</b> = vec(<b>X</b>)}}{\eqn{\textbf{x}=vec(\textbf{X})}{x = vec(X)}}, formed by stacking columns of
\ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}}, is a vector of length 
\ifelse{html}{\out{<i>n</i><i>p</i>}}{\eqn{n p}{n p}}
that also follow a Gaussian distribution with mean 
\ifelse{html}{\out{vec(<b>M</b>)}}{\eqn{vec(\textbf{M})}{vec(M)}} and (co)variance covariance matrix taking the Kronecker form

\ifelse{html}{\out{<p style='text-align:center'><b>&Sigma;</b>&otimes;<b>K</b></p>}}{\deqn{\boldsymbol{\Sigma}\otimes\textbf{K}}{Kronecker(Sigma,K)}}

In the uni-variate case, the problem of near-singularity can be alleviated by penalizing the variance matrix
\ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} by adding positive elements 
\ifelse{html}{\out{&theta;}}{\eqn{\theta}{theta}} to its diagonal, i.e., 
\ifelse{html}{\out{<b>K</b> + &theta;<b>I</b>}}{\eqn{\textbf{K}+\theta\textbf{I}}{K + theta*I}}, where 
\ifelse{html}{\out{<b>I</b>}}{\eqn{\textbf{I}}{I}} is an identity matrix. The same can be applied to the multi-variate case where the Kronecker product (co)variance matrix is penalized with 
\ifelse{html}{\out{<b>&Theta;</b>=\{&theta;<sub>ij</sub>\}}}{\eqn{\boldsymbol{\Theta}=\{\theta_{ij}\}}{Theta=\{theta_ij\}}} of dimensions \ifelse{html}{\out{<i>p</i> &times; <i>p</i>}}{\eqn{p\times p}{p x p}}, where diagonal entries will penalize within feature
\ifelse{html}{\out{<i>i</i>}}{\eqn{i}{i}} and off-diagonals will penalize between features \ifelse{html}{\out{<i>i</i>}}{\eqn{i}{i}}
and \ifelse{html}{\out{<i>j</i>}}{\eqn{j}{j}}. This is,

\ifelse{html}{\out{<p style='text-align:center'><b>&Sigma;</b>&otimes;<b>K</b> + <b>&Theta;</b>&otimes;<b>I</b></p>}}{\deqn{\boldsymbol{\Sigma}\otimes\textbf{K}+\boldsymbol{\Theta}\otimes\textbf{I}}{Kronecker(Sigma,K) + Kronecker(Theta,I)}}

The second Kronecker summand 
\ifelse{html}{\out{<b>&Theta;</b>&otimes;<b>I</b>}}{\eqn{\boldsymbol{\Theta}\otimes\textbf{I}}{Kronecker(Theta,I)}}
is a sparse matrix consisting of non-zero diagonal and sub-diagonals. The \code{Kronecker_cov} function derives the penalized Kronecker (co)variance matrix by computing densely only the first Kronecker summand
\ifelse{html}{\out{<b>&Sigma;</b>&otimes;<b>K</b>}}{\eqn{\boldsymbol{\Sigma}\otimes\textbf{K}}{Kronecker(Sigma,K)}},
and then calculating and adding accordingly only the non-zero entries of 
\ifelse{html}{\out{<b>&Theta;</b>&otimes;<b>I</b>}}{\eqn{\boldsymbol{\Theta}\otimes\textbf{I}}{Kronecker(Theta,I)}}.

\emph{Note}: Swapping the order of the matrices in the above Kronecker operations will yield a different result. In this case the penalized matrix 

\ifelse{html}{\out{<p style='text-align:center'><b>K</b>&otimes;<b>&Sigma;</b> + <b>I</b>&otimes;<b>&Theta;</b></p>}}{\deqn{\textbf{K}\otimes\boldsymbol{\Sigma}+\textbf{I}\otimes\boldsymbol{\Theta}}{Kronecker(K,Sigma) + Kronecker(I,Theta)}}

corresponds to the penalized multi-variate (co)variance matrix of the transposed of the above multi-variate random matrix 
\ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}}, now with features in rows and subjects in columns. This can be achieved by setting \code{swap=TRUE} in the \code{Kronecker_cov} function.

\ifelse{html}{\out{<b>Hadamard product form.</b>}}{\eqn{\textbf{Hadamard product form.}}{Hadamard product form.}}

Assume the random variable \ifelse{html}{\out{<b>x</b><sub>0</sub>}}{\eqn{\textbf{x}_{0}}{x0}} is a subset of \ifelse{html}{\out{<b>x</b>}}{\eqn{\textbf{x}}{x}} containing entries corresponding to specific combinations of subjects and features, then the (co)variance matrix of the vector \ifelse{html}{\out{<b>x</b><sub>0</sub>}}{\eqn{\textbf{x}_{0}}{x0}} will be a Hadamard product formed by the entry-wise product of only the elements of 
\ifelse{html}{\out{<b>&Sigma;</b>}}{\eqn{\boldsymbol{\Sigma}}{Sigma}} and
\ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} involved in the combinations contained in
\ifelse{html}{\out{<b>x</b><sub>0</sub>}}{\eqn{\textbf{x}_{0}}{x0}}; this is

\ifelse{html}{\out{<p style='text-align:center'>(<b>Z</b><sub>1</sub> <b>&Sigma;</b> <b>Z'</b><sub>1</sub>) &odot; (<b>Z</b><sub>2</sub> <b>K</b> <b>Z'</b><sub>2</sub>)</p>}}{\deqn{(\textbf{Z}_{1}\boldsymbol{\Sigma}\textbf{Z}'_{1})\odot(\textbf{Z}_{2}\textbf{K}\textbf{Z}'_{2})}{Hadamard(Z1 Sigma Z1',Z2 K Z2')}}

where
\ifelse{html}{\out{<b>Z</b><sub>1</sub>}}{\eqn{\textbf{Z}_{1}}{Z1}} and
\ifelse{html}{\out{<b>Z</b><sub>2</sub>}}{\eqn{\textbf{Z}_{2}}{Z2}} are incidence matrices mapping from entries of the random variable \ifelse{html}{\out{<b>x</b><sub>0</sub>}}{\eqn{\textbf{x}_{0}}{x0}} to rows (and columns) of \ifelse{html}{\out{<b>&Sigma;</b>}}{\eqn{\boldsymbol{\Sigma}}{Sigma}} and \ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}}, respectively. This (co)variance matrix can be obtained using matrix indexing (see \code{help(Hadamard)}), as 

\ifelse{html}{\out{<p style='text-align:center;font-family:courier'>Sigma[IDS,IDS]*K[IDK,IDK]</p>}}{\deqn{\code{Sigma[IDS,IDS]*K[IDK,IDK]}}{Sigma[IDS,IDS]*K[IDK,IDK]}}

where \code{IDS} and \code{IDK} are integer vectors whose entries are the row (and column) number of 
\ifelse{html}{\out{<b>&Sigma;</b>}}{\eqn{\boldsymbol{\Sigma}}{Sigma}} and \ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}}, respectively, that are mapped at each row of
\ifelse{html}{\out{<b>Z</b><sub>1</sub>}}{\eqn{\textbf{Z}_{1}}{Z1}} and
\ifelse{html}{\out{<b>Z</b><sub>2</sub>}}{\eqn{\textbf{Z}_{2}}{Z2}}, respectively.

The penalized version of this Hadamard product (co)variance matrix will be

\ifelse{html}{\out{<p style='text-align:center'>(<b>Z</b><sub>1</sub> <b>&Sigma;</b> <b>Z'</b><sub>1</sub>) &odot; (<b>Z</b><sub>2</sub> <b>K</b> <b>Z'</b><sub>2</sub>) + (<b>Z</b><sub>1</sub> <b>&Theta;</b> <b>Z'</b><sub>1</sub>) &odot; (<b>Z</b><sub>2</sub> <b>I</b> <b>Z'</b><sub>2</sub>)</p>}}{\deqn{(\textbf{Z}_{1}\boldsymbol{\Sigma}\textbf{Z}'_{1})\odot(\textbf{Z}_{2}\textbf{K}\textbf{Z}'_{2})+(\textbf{Z}_{1}\boldsymbol{\Theta}\textbf{Z}'_{1})\odot(\textbf{Z}_{2}\textbf{I}\textbf{Z}'_{2})}{Hadamard(Z1 Sigma Z1',Z2 K Z2') + Hadamard(Z1 Theta Z1',Z2 I Z2')}}

The \code{Hadamard_cov} function derives this penalized (co)variance matrix using matrix indexing, as 

\ifelse{html}{\out{<p style='text-align:center;font-family:courier'>Sigma[IDS,IDS]*K[IDK,IDK] + Theta[IDS,IDS]*I[IDK,IDK]</p>}}{\deqn{\code{Sigma[IDS,IDS]*K[IDK,IDK] + Theta[IDS,IDS]*I[IDK,IDK]}}{Sigma[IDS,IDS]*K[IDK,IDK] + Theta[IDS,IDS]*I[IDK,IDK]}}

Likewise, this function computes densely only the first Hadamard summand
and then calculates and adds accordingly only the non-zero entries of the second summand.

}
\examples{
  require(tensorEVD)
  
  # Random matrix witn n subjects in rows and p features in columns
  n = 20;  p = 5
  X = matrix(rnorm(n*p), ncol=p)
  
  # Variance matrix among rows/columns
  K = tcrossprod(X)      # for rows
  Sigma = crossprod(X)   # for columns
  Theta = crossprod(matrix(rnorm(p*p), ncol=p))
  dim(K)      # n x n matrix
  dim(Sigma)  # p x p matrix
  dim(Theta)  # p x p matrix
  
  # ==============================================
  # Kronecker covariance
  # ==============================================
  G = Kronecker_cov(Sigma, K, Theta = Theta)

  # it must equal to:
  I0 = diag(n)    # diagonal matrix of dimension n
  G0 = Kronecker(Sigma, K) + Kronecker(Theta, I0)
  all.equal(G,G0)
  
  # Swapping the order of the matrices
  G = Kronecker_cov(Sigma, K, Theta, swap = TRUE)

  # in this case the kronecker is swapped:
  G0 = Kronecker(K, Sigma) + Kronecker(I0, Theta)
  all.equal(G,G0)
  \donttest{
  # ==============================================
  # Extra: Selecting specific entries of the output
  n = 150;  p = 120
  X = matrix(rnorm(n*p), ncol=p)
  K = tcrossprod(X)      
  Sigma = crossprod(X)   
  Theta = crossprod(matrix(rnorm(p*p), ncol=p))
  
  # We want only some rows and columns
  rows = c(1,3,5)
  cols = c(10,30,50)
  G = Kronecker_cov(Sigma, K, Theta, rows=rows, cols=cols)

  # this is preferable instead of:
  # G0 = (Kronecker(Sigma, K) + Kronecker(Theta, diag(n)))[rows,cols]
  # all.equal(G,G0)
  }
  # ==============================================
  # Hadamard covariance
  # ==============================================
  IDS = sample(1:p, 500, replace=TRUE)
  IDK = sample(1:n, 500, replace=TRUE)
  
  G = Hadamard_cov(Sigma, K, Theta, IDS=IDS, IDK=IDK)
  
  # it must equal to:
  I0 = diag(n)    # diagonal matrix of dimension n
  G0 = Sigma[IDS,IDS]*K[IDK,IDK] + Theta[IDS,IDS]*I0[IDK,IDK]
  all.equal(G,G0)
}

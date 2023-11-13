
\name{Kronecker covariance}
\alias{Kronecker_cov}
\title{Kronecker variance matrix penalization}
\usage{
Kronecker_cov(K, Sigma = NULL, Theta, byrow = FALSE,
              rows = NULL, cols = NULL, drop = TRUE)
              
}
\arguments{
\item{K}{(numeric) Variance matrix among subjects}

\item{Sigma}{(numeric) A variance matrix among features. Default \code{Sigma=NULL} will consider an identity matrix with the same dimension as \code{Theta}}

\item{Theta}{(numeric) A diagonal-shifting parameter, value to be added to the diagonals of the Kronecker variance matrix. It can be a (symmetric) matrix with the same dimension as \code{Sigma} for within (diagonal) and between (off-diagonal) features shifting}

\item{byrow}{(logical) If \code{FALSE} (default) the output Kronecker covariance matrix corresponds to a vectorized random matrix stacked by columns, otherwise, it is assumed to be stacked by rows}

\item{rows}{(integer) Index which rows of the Kronecker variance are to be returned. Default \code{rows=NULL} will return all the rows}

\item{cols}{(integer) Index which columns of the Kronecker variance are to be returned. Default \code{cols=NULL} return all the columns}

\item{drop}{Either \code{TRUE} or \code{FALSE} to whether return a uni-dimensional vector when output is a matrix with either 1 row or 1 column as per the \code{rows} and \code{cols} arguments}

}
\value{
Returns the penalized Kronecker covariance matrix. It can be a sub-matrix of it as per the \code{rows} and \code{cols} arguments.
}
\description{
Ridge penalization of a Kronecker covariance matrix
}
\details{
Assume that a multi-variate random matrix \ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}} with 
\ifelse{html}{\out{<i>n</i>}}{\eqn{n}{n}} subjects in rows and \ifelse{html}{\out{<i>p</i>}}{\eqn{p}{p}} features in columns follows a matrix Gaussian distribution with certain matrix of means 
\ifelse{html}{\out{<b>M</b>}}{\eqn{\textbf{M}}{M}} and variance-covariance matrix
\ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} of dimension 
\ifelse{html}{\out{<i>n</i> &times; <i>n</i>}}{\eqn{n\times n}{n x n}} between subjects, and 
\ifelse{html}{\out{<b>&Sigma;</b>}}{\eqn{\boldsymbol{\Sigma}}{Sigma}} of dimension 
\ifelse{html}{\out{<i>p</i> &times; <i>p</i>}}{\eqn{p\times p}{p x p}} between features, then its vectorized form 
\ifelse{html}{\out{vec(<b>X</b>)}}{\eqn{vec(\textbf{X})}{vec(X)}}
will also follow a Gaussian distribution with mean 
\ifelse{html}{\out{vec(<b>M</b>)}}{\eqn{vec(\textbf{M})}{vec(M)}} and variance covariance matrix equal to the Kronecker

\ifelse{html}{\out{<p style='text-align:center'><b>&Sigma;</b>&otimes;<b>K</b></p>}}{\deqn{\boldsymbol{\Sigma}\otimes\textbf{K}}{Kronecker(Sigma,K)}}

if the random matrix is vectorized column-wise or 

\ifelse{html}{\out{<p style='text-align:center'><b>K</b>&otimes;<b>&Sigma;</b></p>}}{\deqn{\textbf{K}\otimes\boldsymbol{\Sigma}}{Kronecker(K,Sigma)}}

if the random matrix is vectorized row-wise.

In the uni-variate case, the problem of near-singularity can be alleviated by penalizing the variance matrix
\ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} by adding positive elements 
\ifelse{html}{\out{&theta;}}{\eqn{\theta}{theta}} to its diagonal, i.e., 
\ifelse{html}{\out{<b>K</b> + &theta;<b>I</b>}}{\eqn{\textbf{K}+\theta\textbf{I}}{K + theta*I}}, where 
\ifelse{html}{\out{<b>I</b>}}{\eqn{\textbf{I}}{I}} is an identity matrix. The same can be applied to the multi-variate case where the Kronecker variance matrix is penalized with 
\ifelse{html}{\out{<b>&Theta;</b>=\{&theta;<sub>ij</sub>\}}}{\eqn{\boldsymbol{\Theta}=\{\theta_{ij}\}}{Theta=\{theta_ij\}}} of dimensions \ifelse{html}{\out{<i>p</i> &times; <i>p</i>}}{\eqn{p\times p}{p x p}}, where diagonal entries will penalize within feature
\ifelse{html}{\out{<i>i</i>}}{\eqn{i}{i}} and off-diagonals will penalize between features \ifelse{html}{\out{<i>i</i>}}{\eqn{i}{i}}
and \ifelse{html}{\out{<i>j</i>}}{\eqn{j}{j}}. This is, 

\ifelse{html}{\out{<p style='text-align:center'><b>&Sigma;</b>&otimes;<b>K</b> + <b>&Theta;</b>&otimes;<b>I</b></p>}}{\deqn{\boldsymbol{\Sigma}\otimes\textbf{K}+\boldsymbol{\Theta}\otimes\textbf{I}}{Kronecker(Sigma,K) + Kronecker(Theta,I)}}

if the random matrix is vectorized column-wise or 

\ifelse{html}{\out{<p style='text-align:center'><b>K</b>&otimes;<b>&Sigma;</b> + <b>I</b>&otimes;<b>&Theta;</b></p>}}{\deqn{\textbf{K}\otimes\boldsymbol{\Sigma}+\textbf{I}\otimes\boldsymbol{\Theta}}{Kronecker(K,Sigma) + Kronecker(I,Theta)}}

if the random matrix is vectorized row-wise.

Specific rows and columns from this Kronecker can be obtained as per the \code{rows} and \code{cols} arguments without forming the whole Kronecker product (see \code{help(Kronecker)}).

}
\examples{
  require(tensorEVD)
  
  # Random matrix witn n subjects in rows and p features in columns
  n = 20
  p = 5
  X = matrix(rnorm(n*p), ncol=p)
  
  # Variance matrix among rows/columns
  K = tcrossprod(X)      # for rows
  Sigma = crossprod(X)   # for columns
  dim(K)      # n x n matrix
  dim(Sigma)  # p x p matrix
  
  # Several examples of penalizing the Kronecker variance
  # ==============================================
  # Example 1. Add unique value 
  # ==============================================
  theta = 1.0
  G = Kronecker_cov(K, Sigma, Theta = theta)

  # it must equal to:
  I0 = diag(n)    # diagonal matrix of dimension n
  Theta0 = matrix(theta, nrow=p, ncol=p)
  G0 = kronecker(Sigma, K) + kronecker(Theta0, I0)
  max(abs(G-G0))
  
  # ==============================================
  # Example 2. Add feature-specific value 
  # ==============================================
  theta = rnorm(p)^2    # One value for each feature
  G = Kronecker_cov(K, Sigma, Theta = theta)

  # it must equal to:
  Theta0 = diag(theta)
  G0 = kronecker(Sigma, K) + kronecker(Theta0, I0)
  max(abs(G-G0))
  
  # ==============================================
  # Example 3. Add specific values within same feature
  #            and between different features 
  # ==============================================
  Theta = crossprod(matrix(rnorm(p*p), ncol=p))
  G = Kronecker_cov(K, Sigma, Theta = Theta)

  # it must equal to:
  G0 = kronecker(Sigma, K) + kronecker(Theta, I0)
  max(abs(G-G0))
  
  # Assume that random matrix X is stacked row-wise
  G = Kronecker_cov(K, Sigma, Theta = Theta, byrow = TRUE)

  # in this case the kronecker is inverted:
  G0 = kronecker(K, Sigma) + kronecker(I0, Theta)
  max(abs(G-G0))
  \donttest{
  # ==============================================
  # Extra: Selecting specific entries of the output
  # ==============================================
  n = 150
  p = 120
  X = matrix(rnorm(n*p), ncol=p)
  K = tcrossprod(X)      
  Sigma = crossprod(X)   
  Theta = crossprod(matrix(rnorm(p*p), ncol=p))
  
  # We want only some rows and columns
  rows = c(1,3,5)
  cols = c(10,30,50)
  G = Kronecker_cov(K, Sigma, Theta = Theta, rows=rows, cols=cols)

  # this is preferable instead of:
  # I0 = diag(n)
  # G0 = (kronecker(Sigma, K) + kronecker(Theta, I0))[rows,cols]
  # max(abs(G-G0))
  }

}


\name{Kronecker product}
\alias{Kronecker}
\title{Kronecker product}
\usage{
Kronecker(A, B, rows = NULL, cols = NULL, 
          make.dimnames = FALSE, drop = TRUE) 

}
\arguments{
\item{A}{(numeric) Left numeric matrix}

\item{B}{(numeric) Right numeric matrix}

\item{rows}{(integer) Index which rows of the Kronecker are to be returned. They must range from 1 to \code{nrow(A)*nrow(B)}. Default \code{rows=NULL} will return all the rows}

\item{cols}{(integer) Index which columns of the Kronecker are to be returned. They must range from 1 to \code{ncol(A)*ncol(B)}. Default \code{cols=NULL} return all the columns}

\item{drop}{Either \code{TRUE} or \code{FALSE} to whether return a uni-dimensional vector when output is a matrix with either 1 row or 1 column as per the \code{rows} and \code{cols} arguments}

\item{make.dimnames}{\code{TRUE} or \code{FALSE} to whether add rownames and colnames attributes to the output}
}
\value{
Returns the Kronecker product matrix. It can be a sub-matrix of it as per the \code{rows} and \code{cols} arguments.
}
\description{
Computes the direct Kronecker product between two matrices
}
\details{
For any two matrices
\ifelse{html}{\out{<b>A</b>=\{a<sub>ij</sub>\}}}{\eqn{\textbf{A}=\{a_{ij}\}}{A=\{a_ij\}}} of dimensions 
\ifelse{html}{\out{<i>m</i> &times; <i>n</i>}}{\eqn{m\times n}{m x n}} and
\ifelse{html}{\out{<b>B</b>=\{b<sub>ij</sub>\}}}{\eqn{\textbf{B}=\{b_{ij}\}}{B=\{b_ij\}}} of dimensions
\ifelse{html}{\out{<i>p</i> &times; <i>q</i>}}{\eqn{p\times q}{p x q}}, 
the direct Kronecker product between them is a matrix defined as the block matrix

\ifelse{html}{\out{<p style='text-align:center'><b>A</b>&otimes;<b>B</b> = \{a<sub>ij</sub><b>B</b>\}</p>}}{\deqn{\textbf{A}\otimes\textbf{B}=\{a_{ij}\textbf{B}\}}{Kronecker(A,B) = \{a_ij B\}}}

which is of dimensions \ifelse{html}{\out{<i>mp</i> &times; <i>nq</i>}}{\eqn{mp\times nq}{mp x nq}}.

Selecting specific rows and columns from the Kronecker can be done by pre- and post- multiplication with incidence matrices

\ifelse{html}{\out{<p style='text-align:center'><b>R</b> (<b>A</b>&otimes;<b>B</b>) <b>C'</b></p>}}{\deqn{\textbf{R}(\textbf{A}\otimes\textbf{B})\textbf{C}'}{R (Kronecker(A,B)) C'}}

where
\ifelse{html}{\out{<b>R</b>}}{\eqn{\textbf{R}}{R}} is an incidence matrix for rows that can be formed by an integer vector \code{rows}, and 
\ifelse{html}{\out{<b>C</b>}}{\eqn{\textbf{C}}{C}} is an incidence matrix for columns that can be formed by an
integer vector \code{cols}.
This sub-matrix of the Kronecker can be obtained by matrix indexing using \code{rows} and \code{cols} as

\ifelse{html}{\out{<p style='text-align:center;font-family:courier'>Kronecker(A,B)[rows,cols]</p>}}{\deqn{\code{Kronecker(A,B)[rows,cols]}}{Kronecker(A,B)[rows,cols]}}

The function computes this sub-matrix of the Kronecker product directly from \ifelse{html}{\out{<b>A</b>}}{\eqn{\textbf{A}}{A}} and \ifelse{html}{\out{<b>B</b>}}{\eqn{\textbf{B}}{B}} without forming the whole Kronecker product. This is very useful if a relatively small number of row/columns are to be selected. 

}
\examples{
  require(tensorEVD)
  
  # Kronecker product of 2 vectors
  A = rnorm(3)
  B = rnorm(4) 
  K = Kronecker(A, B)
  # it must equal when using from the R-base package:
  K0 = kronecker(A, B)
  max(abs(K-K0))
  
  # Kronecker product of 2 matrices
  A = matrix(rnorm(12), ncol=3)
  B = matrix(rnorm(4), ncol=2)
  K = Kronecker(A, B)
  # (it must equal (but faster) to:)
  K0 = kronecker(A, B)
  max(abs(K-K0))
  \donttest{
  # Subsetting rows/columns from the Kronecker
  A = matrix(rnorm(100*150), ncol=150)
  B = matrix(rnorm(100*120), ncol=120)
  rows = c(1,3,5,7)
  cols = c(10,20,30,50)
  K = Kronecker(A, B, rows=rows, cols=cols)
  # (it must equal (but much faster) to:)
  K0 = kronecker(A, B)[rows,cols]
  max(abs(K-K0))
  }

}

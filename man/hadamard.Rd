
\name{Hadamard product}
\alias{Hadamard}
\title{Hadamard product}
\usage{
Hadamard(A, B, rowsA, rowsB,
         colsA = NULL, colsB = NULL,
         make.dimnames = FALSE,
         drop = TRUE)

}
\arguments{
\item{A}{(numeric) Left numeric matrix}

\item{B}{(numeric) Right numeric matrix}

\item{rowsA}{(character/integer) Vector of length \ifelse{html}{\out{<i>m</i>}}{\eqn{m}{m}} with either names or indices mapping from rows of \code{A} into the resulting hadamard product}

\item{rowsB}{(character/integer) Vector of length \ifelse{html}{\out{<i>m</i>}}{\eqn{m}{m}} with names or indices mapping from rows of \code{B} into the resulting hadamard product}

\item{colsA}{(character/integer) (Optional) Similar to \code{rowsA}, vector of length \ifelse{html}{\out{<i>n</i>}}{\eqn{n}{n}} for columns. If \code{colsA = NULL}, it is assumed to be equal to \code{rowsA}}

\item{colsB}{(character/integer) (Optional) Similar to \code{rowsB}, vector of length \ifelse{html}{\out{<i>n</i>}}{\eqn{n}{n}} for columns. If \code{colsB = NULL}, it is assumed to be equal to \code{rowsB}}

\item{drop}{Either \code{TRUE} or \code{FALSE} to whether return a uni-dimensional vector when output is a matrix with either 1 row or 1 column as per the \code{rows} and \code{cols} arguments}

\item{make.dimnames}{\code{TRUE} or \code{FALSE} to whether add rownames and colnames attributes to the output}
}
\value{
Returns a matrix containing the Hadamard product.
}
\description{
Computes the Hadamard product between two matrices
}
\details{
Computes the \ifelse{html}{\out{<i>m</i> &times; <i>n</i>}}{\eqn{m\times n}{m x n}} Hadamard product (aka element-wise or entry-wise product) matrix between matrices 
\ifelse{html}{\out{<b>A</b><sub>0</sub> = <b>R</b><sub>1</sub> <b>A</b> <b>C'</b><sub>1</sub>}}{\eqn{\textbf{A}_{0}=\textbf{R}_{1}\textbf{A}\textbf{C}'_{1}}{A0 = R1 A C1'}} and
\ifelse{html}{\out{<b>B</b><sub>0</sub> = <b>R</b><sub>2</sub> <b>B</b> <b>C'</b><sub>2</sub>}}{\eqn{\textbf{B}_{0}=\textbf{R}_{2}\textbf{B}\textbf{C}'_{2}}{B0 = R2 B C2'}},

\ifelse{html}{\out{<p style='text-align:center'>(<b>R</b><sub>1</sub> <b>A</b> <b>C'</b><sub>1</sub>) &odot; (<b>R</b><sub>2</sub> <b>B</b> <b>C'</b><sub>2</sub>)</p>}}{\deqn{(\textbf{R}_{1}\textbf{A}\textbf{C}'_{1})\odot(\textbf{R}_{2}\textbf{B}\textbf{C}'_{2})}{Hadamard(R1 A C1', R2 B C2')}}

where
\ifelse{html}{\out{<b>R</b><sub>1</sub>}}{\eqn{\textbf{R}_{1}}{R1}} and
\ifelse{html}{\out{<b>R</b><sub>2</sub>}}{\eqn{\textbf{R}_{2}}{R2}} are incidence matrices for rows that can be formed by integer vectors
\code{rowsA} and \code{rowsB} of length \ifelse{html}{\out{<i>m</i>}}{\eqn{m}{m}}, respectively, and 
\ifelse{html}{\out{<b>C</b><sub>1</sub>}}{\eqn{\textbf{C}_{1}}{C1}} and
\ifelse{html}{\out{<b>C</b><sub>2</sub>}}{\eqn{\textbf{C}_{2}}{C2}} are incidence matrices for columns that can be formed by integer vectors
\code{colsA} and \code{colsB} of length \ifelse{html}{\out{<i>n</i>}}{\eqn{n}{n}}, respectively.

Matrices \ifelse{html}{\out{<b>A</b><sub>0</sub>}}{\eqn{\textbf{A}_{0}}{A0}} and 
\ifelse{html}{\out{<b>B</b><sub>0</sub>}}{\eqn{\textbf{B}_{0}}{B0}} can be obtained by matrix indexing as \code{A[rowsA,colsA]} and \code{B[rowsB,colsB]}, respectively. Therefore, the Hadamard product can be obtained directly as

\ifelse{html}{\out{<p style='text-align:center;font-family:courier'>A[rowsA,colsA]*B[rowsB,colsB]</p>}}{\deqn{\code{A[rowsA,colsA]*B[rowsB,colsB]}}{A[rowsA,colsA]*B[rowsB,colsB]}}

The function computes the Hadamard product directly from \ifelse{html}{\out{<b>A</b>}}{\eqn{\textbf{A}}{A}} and \ifelse{html}{\out{<b>B</b>}}{\eqn{\textbf{B}}{B}} without forming \ifelse{html}{\out{<b>A</b><sub>0</sub>}}{\eqn{\textbf{A}_{0}}{A0}} or \ifelse{html}{\out{<b>B</b><sub>0</sub>}}{\eqn{\textbf{B}_{0}}{B0}} matrices.

}
\examples{
  require(tensorEVD)
  
  # ==============================================
  # Example 1. Indexing using integers
  # ==============================================
  # Generate rectangular matrices A (nrowA x ncolA) and B (nrowB x ncolB)
  nA = c(10,15)
  nB = c(12,8)
  A = matrix(rnorm(nA[1]*nA[2]), nrow=nA[1])
  B = matrix(rnorm(nB[1]*nB[2]), nrow=nB[1])
  
  # Define size of the Hadamard n1 x n2
  n1 = 1000
  n2 = 500
  rowsA = sample(seq(nA[1]), n1, replace=TRUE)
  rowsB = sample(seq(nB[1]), n1, replace=TRUE)
  colsA = sample(seq(nA[2]), n2, replace=TRUE)
  colsB = sample(seq(nB[2]), n2, replace=TRUE)
  
  # Direct hadamard product
  K1 = A[rowsA,colsA]*B[rowsB,colsB]
  
  # Using 'Hadamard' function
  K2 = Hadamard(A, B, rowsA, rowsB, colsA, colsB)
  
  max(K1-K2)  # They should be equal
  
  # ==============================================
  # Example 2. Indexing using row/column names
  # ==============================================
  # Generate squared symmetric matrices A and B 
  nA = 20
  nB = 15
  A = tcrossprod(matrix(rnorm(nA*nA), nrow=nA, dimnames=list(paste0("id",seq(nA)))))
  B = tcrossprod(matrix(rnorm(nB*nB), nrow=nB, dimnames=list(paste0("id",seq(nB)))))
  
  # Define size of the Hadamard n x n
  n = 1000
  IDA = sample(rownames(A), n, replace=TRUE)
  IDB = sample(rownames(B), n, replace=TRUE)
  
  # Direct hadamard product
  K1 = A[IDA,IDA]*B[IDB,IDB]
  
  # Using 'Hadamard' function
  K2 = Hadamard(A, B, IDA, IDB)
  
  max(K1-K2)  # They should be equal

}

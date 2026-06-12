
\name{Tensor EVD}
\alias{tensorEVD}
\title{Tensor EVD}
\usage{
tensorEVD(K1, K2, ID1, ID2, EVD1 = NULL,
          EVD2 = NULL, nu = NULL, alpha = 1,
          d.min = .Machine$double.eps*10,
          norm = NULL, order = NULL, 
          make.dimnames = FALSE, verbose = FALSE)
}
\arguments{
\item{K1, K2}{(numeric) Covariance structure matrices}

\item{ID1}{(character/integer) Vector of length \ifelse{html}{\out{<i>n</i>}}{\eqn{n}{n}} with either names or indices mapping from rows/columns of \code{K1} into the resulting tensor product}

\item{ID2}{(character/integer) Vector of length \ifelse{html}{\out{<i>n</i>}}{\eqn{n}{n}} with either names or indices mapping from rows/columns of \code{K2} into the resulting tensor product}

\item{EVD1}{(list) (Optional) Eigenvalues and eigenvectors of \code{K1} as produced by the \code{eigen} function. If provided, \code{K1} is not needed}

\item{EVD2}{(list) (Optional) Eigenvalues and eigenvectors of \code{K2} as produced by the \code{eigen} function. If provided, \code{K2} is not needed}

\item{nu}{(integer) Number of tensor eigenvectors to be returned}

\item{alpha}{(numeric) Proportion of variance explained by the tensor eigenvectors. Only eigenvectors jointly explaining such proportion of variance are returned}

\item{d.min}{(numeric) Tensor eigenvalue threshold. Default is a numeric zero. Only eigenvectors with eigenvalue passing this threshold are returned}

\item{norm}{(numeric) Vector with (precalculated) tensor eigenvectors' norms. This is used to scale tensor eigenvectors and eigenvalues before sorting. The length of this should equal to \code{length(d1)} times \code{length(d2)}. If \code{NULL}, weights equal to 1 are assumed}

\item{order}{(integer) Vector with indices to sort the eigenvalues product \code{Kronecker(d1,d2)}. The length of this should equal to \code{length(d1)} times \code{length(d2)}. If \code{NULL}, order is internally calculated}

\item{make.dimnames}{\code{TRUE} or \code{FALSE} to whether add \code{rownames} and \code{colnames} attributes to the output}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether show progress}
}
\value{
Returns a list object that contains the elements:
\itemize{
  \item \code{values}: (vector) resulting tensor eigenvalues.
  \item \code{vectors}: (matrix) resulting normalized tensor eigenvectors.
  \item \code{norm}: (vector) norm of the unsorted tensor eigenvectors.
  \item \code{order}: (vector) indices used to sort tensor eigenvectors.
  \item \code{totalVar}: (numeric) total variance of the tensor matrix product.
}
}
\description{
Fast eigen value decomposition (EVD) of the Hadamard product of two matrices
}
\details{
Let the \ifelse{html}{\out{<i>n</i> &times; <i>n</i>}}{\eqn{n\times n}{n x n}} matrix \ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} to be the Hadamard product (aka element-wise or entry-wise product) involving two smaller matrices 
\ifelse{html}{\out{<b>K</b><sub>1</sub>}}{\eqn{\textbf{K}_{1}}{K1}} and
\ifelse{html}{\out{<b>K</b><sub>2</sub>}}{\eqn{\textbf{K}_{2}}{K2}} of dimensions 
\ifelse{html}{\out{<i>n</i><sub>1</sub>}}{\eqn{n_{1}}{n1}} and \ifelse{html}{\out{<i>n</i><sub>2</sub>}}{\eqn{n_{2}}{n2}}, respectively,

\ifelse{html}{\out{<p style='text-align:center'><b>K</b> = (<b>Z</b><sub>1</sub> <b>K</b><sub>1</sub> <b>Z'</b><sub>1</sub>) &odot; (<b>Z</b><sub>2</sub> <b>K</b><sub>2</sub> <b>Z'</b><sub>2</sub>)</p>}}{\deqn{\textbf{K}=(\textbf{Z}_{1}\textbf{K}_{1}\textbf{Z}'_{1})\odot(\textbf{Z}_{2}\textbf{K}_{2}\textbf{Z}'_{2})}{K = Hadamard(Z1 K1 Z1', Z2 K2 Z2')}}

where
\ifelse{html}{\out{<b>Z</b><sub>1</sub>}}{\eqn{\textbf{Z}_{1}}{Z1}} and
\ifelse{html}{\out{<b>Z</b><sub>2</sub>}}{\eqn{\textbf{Z}_{2}}{Z2}} are incidence matrices mapping from rows (and columns) of the resulting Hadamard to rows (and columns) of \ifelse{html}{\out{<b>K</b><sub>1</sub>}}{\eqn{\textbf{K}_{1}}{K1}} and
\ifelse{html}{\out{<b>K</b><sub>2</sub>}}{\eqn{\textbf{K}_{2}}{K2}}, respectively. 

Let the eigenvalue decomposition (EVD) of \ifelse{html}{\out{<b>K</b><sub>1</sub>}}{\eqn{\textbf{K}_{1}}{K1}} and
\ifelse{html}{\out{<b>K</b><sub>2</sub>}}{\eqn{\textbf{K}_{2}}{K2}} to be
\ifelse{html}{\out{<b>K</b><sub>1</sub> = <b>U</b><sub>1</sub> <b>D</b><sub>1</sub> <b>U'</b><sub>1</sub>}}{\eqn{\textbf{K}_{1}=\textbf{U}_{1}\textbf{D}_{1}\textbf{U}'_{1}}{K1 = U1 D1 U1'}} and 
\ifelse{html}{\out{<b>K</b><sub>2</sub> = <b>U</b><sub>2</sub> <b>D</b><sub>2</sub> <b>U'</b><sub>2</sub>}}{\eqn{\textbf{K}_{2}=\textbf{U}_{2}\textbf{D}_{2}\textbf{U}'_{2}}{K2 = U2 D2 U2'}}. 
Using properties of the Hadamard and Kronecker products, an EVD of the Hadamard product
\ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} can be approximated using the EVD of 
\ifelse{html}{\out{<b>K</b><sub>1</sub>}}{\eqn{\textbf{K}_{1}}{K1}} and 
\ifelse{html}{\out{<b>K</b><sub>2</sub>}}{\eqn{\textbf{K}_{2}}{K2}} as

\ifelse{html}{\out{<p style='text-align:center'><b>K = U D U'</b></p>}}{\deqn{\textbf{K}=\textbf{U}\textbf{D}\textbf{U}'}{K = U D U'}}

where \ifelse{html}{\out{<b>D</b> = <b>D</b><sub>1</sub>&otimes;<b>D</b><sub>2</sub>}}{\eqn{\textbf{D}=\textbf{D}_1\otimes\textbf{D}_2}{D = Kronecker(D1,D2)}} is a diagonal matrix containing
\ifelse{html}{\out{<i>N</i> = <i>n</i><sub>1</sub> &times; <i>n</i><sub>2</sub>}}{\eqn{N=n_{1}\times n_{2}}{N = n1 x n2}} tensor eigenvalues 
\ifelse{html}{\out{d<sub>1</sub> &ge; ... &ge; d<sub>N</sub> &ge; 0}}{\eqn{d_1\ge...\ge d_N \ge 0}{d_1 >= ...>= d_N >= 0}}, and
\ifelse{html}{\out{<b>U</b> = (<b>Z</b><sub>1</sub>&Star;<b>Z</b><sub>2</sub>)(<b>U</b><sub>1</sub>&otimes;<b>U</b><sub>2</sub>) = [<b>u</b><sub>1</sub>,...,<b>u</b><sub>N</sub>]}}{\eqn{\textbf{U}=(\textbf{Z}_1\star\textbf{Z}_2)(\textbf{U}_1\otimes\textbf{U}_2)=[\textbf{u}_1,...,\textbf{u}_N]}{U = (Z1*Z2)(Kronecker(U1,U2)) = [u_1,...,u_N]}} is matrix containing \ifelse{html}{\out{<i>N</i>}}{\eqn{N}{N}} tensor eigenvectors
\ifelse{html}{\out{<b>u</b><sub>k</sub>}}{\eqn{\textbf{u}_{k}}{u_k}}. 
The term 
\ifelse{html}{\out{<b>Z</b><sub>1</sub>&Star;<b>Z</b><sub>2</sub>}}{\eqn{\textbf{Z}_1\star\textbf{Z}_2}{Z1*Z2}} is the 
"face-splitting product" (aka "transposed Khatri–Rao product") of matrices 
\ifelse{html}{\out{<b>Z</b><sub>1</sub>}}{\eqn{\textbf{Z}_{1}}{Z1}} and
\ifelse{html}{\out{<b>Z</b><sub>2</sub>}}{\eqn{\textbf{Z}_{2}}{Z2}}, obtained as a row-by-row Kronecker product.

Each tensor eigenvector \ifelse{html}{\out{<i>k</i>}}{\eqn{k}{k}} is derived separately as a Hadamard product using the corresponding 
\ifelse{html}{\out{<i>i(k)</i>}}{\eqn{i(k)}{i(k)}} and \ifelse{html}{\out{<i>j(k)</i>}}{\eqn{j(k)}{j(k)}} eigenvectors 
\ifelse{html}{\out{<b>u</b><sub>1i(k)</sub>}}{\eqn{\textbf{u}_{1i(k)}}{u_1i(k)}} and
\ifelse{html}{\out{<b>u</b><sub>2j(k)</sub>}}{\eqn{\textbf{u}_{2j(k)}}{u_2j(k)}} of 
\ifelse{html}{\out{<b>U</b><sub>1</sub>}}{\eqn{\textbf{U}_{1}}{U1}} and 
\ifelse{html}{\out{<b>U</b><sub>2</sub>}}{\eqn{\textbf{U}_{2}}{U2}}, respectively, that form the \ifelse{html}{\out{<i>k</i>}}{\eqn{k}{k}} eigenvector of 
\ifelse{html}{\out{<b>U</b><sub>1</sub>&otimes;<b>U</b><sub>2</sub>}}{\eqn{\textbf{U}_1\otimes\textbf{U}_2}{Kronecker(U1,U2)}}, this is

\ifelse{html}{\out{<p style='text-align:center'><b>u</b><sub>k</sub> = (<b>Z</b><sub>1</sub><b>u</b><sub>1i(k)</sub>)&odot;(<b>Z</b><sub>2</sub><b>u</b><sub>2j(k)</sub>)</p>}}{\deqn{\textbf{u}_{k}=(\textbf{Z}_{1}\textbf{u}_{1i(k)})\odot(\textbf{Z}_{2}\textbf{u}_{2j(k)})}{v_k = Hadamard(Z1 u_1i(k), Z2 u_2j(k))}}

The \code{tensorEVD} function derives each of these eigenvectors \ifelse{html}{\out{<b>u</b><sub>k</sub>}}{\eqn{\textbf{u}_{k}}{u_k}} by matrix indexing using integer vectors \code{ID1} and \code{ID2}. The entries of these vectors are the row (and column) number of 
\ifelse{html}{\out{<b>K</b><sub>1</sub>}}{\eqn{\textbf{K}_{1}}{K1}} and 
\ifelse{html}{\out{<b>K</b><sub>2</sub>}}{\eqn{\textbf{K}_{2}}{K2}} that are mapped at each row of \ifelse{html}{\out{<b>Z</b><sub>1</sub>}}{\eqn{\textbf{Z}_{1}}{Z1}} and
\ifelse{html}{\out{<b>Z</b><sub>2</sub>}}{\eqn{\textbf{Z}_{2}}{Z2}}, respectively.

}
\examples{
  require(tensorEVD)
  set.seed(195021)
  
  # Generate matrices K1 and K2 of dimensions n1 and n2
  n1 = 10; n2 = 15
  K1 = crossprod(matrix(rnorm(n1*(n1+10)), ncol=n1))
  K2 = crossprod(matrix(rnorm(n2*(n2+10)), ncol=n2))

  # (a) Example 1. Full design (Kronecker product)
  ID1 = rep(seq(n1), each=n2)
  ID2 = rep(seq(n2), times=n1)

  # Direct EVD of the Hadamard product
  K = K1[ID1,ID1]*K2[ID2,ID2]
  EVD0 = eigen(K)

  # Tensor EVD using K1 and K2
  EVD = tensorEVD(K1, K2, ID1, ID2)

  # Eigenvectors and eigenvalues are numerically equal
  all.equal(EVD0$values, EVD$values)
  all.equal(abs(EVD0$vectors), abs(EVD$vectors)) 

  # (b) If a proportion of variance explained is specified, 
  # only the eigenvectors needed to explain such proportion are derived
  alpha = 0.95
  EVD = tensorEVD(K1, K2, ID1, ID2, alpha=alpha)
  dim(EVD$vectors)

  # For the direct EVD
  varexp = cumsum(EVD0$values/sum(EVD0$values))
  index = 1:which.min(abs(varexp-alpha))
  dim(EVD0$vectors[,index])

  # (c) Example 2. Incomplete design (Hadamard product)
  # Eigenvectors and eigenvalues are no longer equivalent
  n = n1*n2   # Sample size n
  ID1 = sample(seq(n1), n, replace=TRUE) # Randomly sample of ID1
  ID2 = sample(seq(n2), n, replace=TRUE) # Randomly sample of ID2

  K = K1[ID1,ID1]*K2[ID2,ID2]
  EVD0 = eigen(K)
  EVD = tensorEVD(K1, K2, ID1, ID2)

  all.equal(EVD0$values, EVD$values)
  all.equal(abs(EVD0$vectors), abs(EVD$vectors)) 
  
  # However, the sum of the eigenvalues is equal to the trace(K)
  c(sum(EVD0$values), sum(EVD$values), sum(diag(K)))

  # And provide the same approximation for K
  K01 = EVD0$vectors\%*\%diag(EVD0$values)\%*\%t(EVD0$vectors)
  K02 = EVD$vectors\%*\%diag(EVD$values)\%*\%t(EVD$vectors)
  c(all.equal(K,K01), all.equal(K,K02))

  # When n is different from N=n1xn2, both methods provide different 
  # number or eigenvectors/eigenvalues. The eigen function provides 
  # a number of eigenvectors equal to the minimum between n and N
  # for the tensorEVD, this number is always N

  # (d) Sample size n being half of n1 x n2
  n = n1*n2/2
  ID1 = sample(seq(n1), n, replace=TRUE)
  ID2 = sample(seq(n2), n, replace=TRUE)

  K = K1[ID1,ID1]*K2[ID2,ID2]
  EVD0 = eigen(K)
  EVD = tensorEVD(K1, K2, ID1, ID2)

  c(eigen=sum(EVD0$values>1E-10), tensorEVD=sum(EVD$values>1E-10))

  # (e) Sample size n being twice n1 x n2
  n = n1*n2*2
  ID1 = sample(seq(n1), n, replace=TRUE)
  ID2 = sample(seq(n2), n, replace=TRUE)

  K = K1[ID1,ID1]*K2[ID2,ID2]
  EVD0 = eigen(K)
  EVD = tensorEVD(K1, K2, ID1, ID2)

  c(eigen=sum(EVD0$values>1E-10), tensorEVD=sum(EVD$values>1E-10))
}


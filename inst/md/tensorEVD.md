### Tensor EVD

Let the <i>n</i> &times; <i>n</i> matrix <b>K</b> to be a Hadamard product involving two smaller matrices <b>K</b><sub>1</sub> and <b>K</b><sub>2</sub> of dimensions <i>n</i><sub>1</sub> and <i>n</i><sub>2</sub>, respectively, 

<p align="center">
<b>K</b> = (<b>Z</b><sub>1</sub> <b>K</b><sub>1</sub> <b>Z'</b><sub>1</sub>) &odot; (<b>Z</b><sub>2</sub> <b>K</b><sub>2</sub> <b>Z'</b><sub>2</sub>)
</p>

where <b>Z</b><sub>1</sub> and <b>Z</b><sub>2</sub> are incidence matrices for <b>K</b><sub>1</sub> and <b>K</b><sub>2</sub>, respectively.

Let the eigenvalue decomposition (EVD) of <b>K</b><sub>1</sub> and <b>K</b><sub>2</sub> to be
<b>K</b><sub>1</sub> = <b>V</b><sub>1</sub> <b>D</b><sub>1</sub> <b>V'</b><sub>1</sub> and 
<b>K</b><sub>2</sub> = <b>V</b><sub>2</sub> <b>D</b><sub>2</sub> <b>V'</b><sub>2</sub>. 
Using properties of the Hadamard and Kronecker products, an EVD of the Hadamard product <b>K</b> can be approximated using the EVD of 
<b>K</b><sub>1</sub> and <b>K</b><sub>2</sub> as

<p align="center">
<b>K = V D V'</b>
</p>

where <b>D</b> = <b>D</b><sub>1</sub>&otimes;<b>D</b><sub>2</sub> is a diagonal matrix containing
<i>N</i> = <i>n</i><sub>1</sub> &times; <i>n</i><sub>2</sub> tensor eigenvalues 
d<sub>1</sub> &ge; ... &ge; d<sub>N</sub> &ge; 0, and
<b>V</b> = (<b>Z</b><sub>1</sub>&Star;<b>Z</b><sub>2</sub>)(<b>V</b><sub>1</sub>&otimes;<b>V</b><sub>2</sub>) = [<b>v</b><sub>1</sub>,...,<b>v</b><sub>N</sub>] is a matrix containing <i>N</i> tensor eigenvectors
<b>v</b><sub>k</sub>; here the term 
<b>Z</b><sub>1</sub>&Star;<b>Z</b><sub>2</sub> is the 
"face-splitting product" (aka "transposed Khatri–Rao product") of matrices 
<b>Z</b><sub>1</sub> and <b>Z</b><sub>2</sub>.

Each tensor eigenvector <b>v</b><sub>k</sub> is derived separately as a Hadamard product using the corresponding 
<i>i(k)</i> and <i>j(k)</i> eigenvectors <b>v</b><sub>1i(k)</sub> and <b>v</b><sub>2j(k)</sub> from 
<b>V</b><sub>1</sub> and <b>V</b><sub>2</sub>, respectively, this is

<p align="center">
<b>v</b><sub>k</sub> = (<b>Z</b><sub>1</sub><b>v</b><sub>1i(k)</sub>)&odot;(<b>Z</b><sub>2</sub><b>v</b><sub>2j(k)</sub>)
</p>

This approach can be implemented with the *tensorEVD*() function using integer vectors `ID1` and `ID2` instead of matrices <b>Z</b><sub>1</sub> and <b>Z</b><sub>2</sub>. The EVD of <b>K</b><sub>1</sub> and <b>K</b><sub>2</sub> can be performed using the *eigen*() R-function.

#### Example 1. Complete design
When <b>Z</b><sub>1</sub> and <b>Z</b><sub>2</sub> are such that the Hadamard <b>K</b> represent a Kronecker product between <b>K</b><sub>1</sub> and <b>K</b><sub>2</sub> (i.e, <i>n</i> = <i>n</i><sub>1</sub> &times; <i>n</i><sub>2</sub> combinations, each element of <b>K</b><sub>1</sub> crossed with one and only one element of <b>K</b><sub>2</sub>), then the results of the *tensorEVD*() are exactly the same as those obtained by performing the EVD directly on <b>K</b>. 
```r
n1 = 20; n2 = 30
K1 = crossprod(matrix(rnorm(n1*(n1+10)), ncol=n1))
K2 = crossprod(matrix(rnorm(n2*(n2+10)), ncol=n2))

ID1 = rep(seq(n1), each=n2)
ID2 = rep(seq(n2), times=n1)

# Direct EVD of the Hadamard product
K = Hadamard(K1, K2, ID1, ID2)   # Same as K = K1[ID1,ID1]*K2[ID2,ID2]
EVD0 = eigen(K)

# Tensor EVD using K1 and K2
EVD1 = eigen(K1)
EVD2 = eigen(K2)
EVD = tensorEVD(V1=EVD1$vectors, d1=EVD1$values, V2=EVD2$vectors, d2=EVD2$values, ID1=ID1, ID2=ID2)

# Eigenvectors and eigenvalues are numerically equal
max(abs(EVD0$values-EVD$values)) # near zero
max(abs(EVD0$vectors)-abs(EVD$vectors)) # near zero
```

Matrices <b>K</b><sub>1</sub> and <b>K</b><sub>2</sub> can be also passed to the function, in this case, the EVD will be computed internally,
```r
EVD = tensorEVD(K1, K2, ID1, ID2)

max(abs(EVD0$values-EVD$values))
max(abs(EVD0$vectors)-abs(EVD$vectors))
```

#### Example 2. Incomplete design
For incomplete (i.e., less than <i>n</i><sub>1</sub> &times; <i>n</i><sub>2</sub> unique combinations) or unbalanced (i.e., combinations appearing at different frequencies) designs, eigenvectors and eigenvalues are no longer equivalent
```r
n = n1*n2   # Sample size n
ID1 = sample(seq(n1), n, replace=TRUE) # Randomly sample of ID1
ID2 = sample(seq(n2), n, replace=TRUE) # Randomly sample of ID2

K = Hadamard(K1, K2, ID1, ID2)
EVD0 = eigen(K)
EVD = tensorEVD(K1, K2, ID1, ID2)

max(abs(EVD0$values-EVD$values))
max(abs(EVD0$vectors)-abs(EVD$vectors))
```

Even the number of eigenvectors is not the same when <i>n</i> &ne; <i>n</i><sub>1</sub> &times; <i>n</i><sub>2</sub>. The *eigen*() function provides an EVD where the number of eigenvectors is equal to the minimum between <i>n</i> and <i>n</i><sub>1</sub> &times; <i>n</i><sub>2</sub>; for the *tensorEVD*() function, this number is always equal to <i>n</i><sub>1</sub> &times; <i>n</i><sub>2</sub>.
```r
# Sample size n being twice n1 x n2
n = n1*n2*2
ID1 = sample(seq(n1), n, replace=TRUE)
ID2 = sample(seq(n2), n, replace=TRUE)

K = Hadamard(K1, K2, ID1, ID2)
EVD0 = eigen(K)
EVD = tensorEVD(K1, K2, ID1, ID2)

c(ncol(EVD0$vectors), ncol(EVD$vectors))

# Sample size n being half of n1 x n2
n = n1*n2/2
ID1 = sample(seq(n1), n, replace=TRUE)
ID2 = sample(seq(n2), n, replace=TRUE)

K = Hadamard(K1, K2, ID1, ID2)
EVD0 = eigen(K)
EVD = tensorEVD(K1, K2, ID1, ID2)

c(ncol(EVD0$vectors), ncol(EVD$vectors))
```
  
However, both methods will produce a total sum of eigenvalues always equal to the <i>trace</i>(<b>K</b>) and will provide the same approximation for 
<b>K = V D V'</b>,
```r
# Sum of eigenvalues
c(sum(EVD0$values), sum(EVD$values), sum(diag(K)))

# Approximation for K
K01 = EVD0$vectors%*%diag(EVD0$values)%*%t(EVD0$vectors)
K02 = EVD$vectors%*%diag(EVD$values)%*%t(EVD$vectors)
c(max(K-K01), max(K-K02))
```

<- [Return](https://github.com/MarcooLopez/tensorEVD/blob/main/README.md)

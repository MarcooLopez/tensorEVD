## ----initialsetup, include=FALSE----------------------------------------------
knitr::opts_chunk$set(cache=FALSE)
library(tensorEVD)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  tensorEVD(K1, K2, ID1, ID2)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Simulating covariance matrices K1 and K2
n1 = 10; n2 = 15
K1 <- crossprod(matrix(rnorm(n1*(n1+10),sd=sqrt(1/n1)), ncol=n1))
K2 <- crossprod(matrix(rnorm(n2*(n2+10),sd=sqrt(1/n2)), ncol=n2))

ID1 <- rep(seq(n1), each=n2)
ID2 <- rep(seq(n2), times=n1)

# Direct EVD of the Hadamard product
K <- Hadamard(K1, K2, ID1, ID2, ID1, ID2)   # Same as K = K1[ID1,ID1]*K2[ID2,ID2]
EVD0 <- eigen(K)

# Tensor EVD using K1 and K2
EVD <- tensorEVD(K1, K2, ID1, ID2)

# Eigenvalues and (absolute) eigenvectors and are numerically equal
all.equal(EVD0$values, EVD$values)
all.equal(abs(EVD0$vectors), abs(EVD$vectors)) 

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
n <- n1*n2   # size of the Hadamard
ID1 <- sample(seq(n1), n, replace=TRUE) # Randomly sample of ID1
ID2 <- sample(seq(n2), n, replace=TRUE) # Randomly sample of ID2

K <- Hadamard(K1, K2, ID1, ID2, ID1, ID2)
EVD0 <- eigen(K)
EVD <- tensorEVD(K1, K2, ID1, ID2)

all.equal(EVD0$values, EVD$values)
all.equal(abs(EVD0$vectors), abs(EVD$vectors)) 

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Sum of eigenvalues
c(sum(EVD0$values), sum(EVD$values), sum(diag(K)))

# Approximation for K
K01 <- EVD0$vectors%*%diag(EVD0$values)%*%t(EVD0$vectors)
K02 <- EVD$vectors%*%diag(EVD$values)%*%t(EVD$vectors)
c(all.equal(K,K01), all.equal(K,K02))

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
n = n1*n2/2    # size of the Hadamard is half of n1 x n2
ID1 <- sample(seq(n1), n, replace=TRUE)
ID2 <- sample(seq(n2), n, replace=TRUE)

K <- Hadamard(K1, K2, ID1, ID2, ID1, ID2)
EVD0 <- eigen(K)
EVD <- tensorEVD(K1, K2, ID1, ID2)

# Number of eigenvectors with positive eigenvalue
c(eigen=sum(EVD0$values>1E-10), tensorEVD=sum(EVD$values>1E-10))

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Size of the Hadamard is three times n1 x n2
# Balanced and replicated case
ID1 <- rep(rep(seq(n1), each=n2), 3)
ID2 <- rep(rep(seq(n2), times=n1), 3)

K <- Hadamard(K1, K2, ID1, ID2, ID1, ID2)
EVD0 <- eigen(K)
EVD <- tensorEVD(K1, K2, ID1, ID2)

c(eigen=sum(EVD0$values>1E-10), tensorEVD=sum(EVD$values>1E-10))

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
alpha <- 0.95
EVD <- tensorEVD(K1, K2, ID1, ID2, alpha=alpha)
ncol(EVD$vectors)

# For the direct EVD
varexp = cumsum(EVD0$values/sum(EVD0$values))
index = 1:which.min(abs(varexp-alpha))
ncol(EVD0$vectors[,index])

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
dimnames(K1) <- list(paste0("i",seq(n1)), paste0("i",seq(n1)))
dimnames(K2) <- list(paste0("j",seq(n2)), paste0("j",seq(n2)))

EVD <- tensorEVD(K1, K2, ID1, ID2, make.dimnames=TRUE)
EVD$vectors[1:6,1:5]

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
EVD2 <- eigen(K2)
EVD0 <- tensorEVD(K1=K1, EVD2=EVD2, ID1=ID1, ID2=ID2)
EVD <- tensorEVD(K1=K1, K2=K2, ID1=ID1, ID2=ID2)

all.equal(EVD0$values, EVD$values)
all.equal(abs(EVD0$vectors), abs(EVD$vectors)) 


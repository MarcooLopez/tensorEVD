## ----initialsetup, include=FALSE----------------------------------------------
knitr::opts_chunk$set(cache=FALSE)
library(tensorEVD)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Simulating matrices A and B
m = 100; n = 150
A <- matrix(rnorm(m*n), ncol=n)
B <- matrix(rnorm(m*n), ncol=n)

# Making the Hadamard product
K1 <- A*B
K2 <- Hadamard(A, B)

all.equal(K1, K2) # should be equal

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
a <- 10
( B <- matrix(1:6, ncol=2) )

a*B    # using the product operator
try(Hadamard(a, B), silent=TRUE)[[1]]
Kronecker(a, B)  # Kronecker instead of Hadamard

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Simulating A and B of different dimensions
m = 100; n = 150
p = 200; q = 120
A <- matrix(rnorm(m*n), ncol=n)  # m x n
B <- matrix(rnorm(p*q), ncol=q)  # p x q

try(A*B, silent=TRUE)[[1]]

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Subsetting rows and columns of A each of length 
# equal to nrow(B) and ncol(B), respectively
rowsA <- sample(seq(nrow(A)), nrow(B), replace=TRUE)
colsA <- sample(seq(ncol(A)), ncol(B), replace=TRUE)

# Making the Hadamard product
K1 <- A[rowsA,colsA]*B
K2 <- Hadamard(A, B, rowsA=rowsA, colsA=colsA)

all.equal(K1, K2)
dim(K2) == dim(B) # has the same dimension as B

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
rowsB <- sample(seq(nrow(B)), nrow(A), replace=TRUE)
colsB <- sample(seq(ncol(B)), ncol(A), replace=TRUE)

K2 <- Hadamard(A, B, rowsB=rowsB, colsB=colsB)
dim(K2) == dim(A) 

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
dm <- c(1000, 2000)  # a Hadamard of 1000 x 2000

# Obtaining a submatrix from A
rowsA <- sample(seq(nrow(A)), dm[1], replace=TRUE)
colsA <- sample(seq(ncol(A)), dm[2], replace=TRUE)

# Obtaining a submatrix from B
rowsB <- sample(seq(nrow(B)), dm[1], replace=TRUE)
colsB <- sample(seq(ncol(B)), dm[2], replace=TRUE)

# Making the Hadamard product
K1 <- A[rowsA,colsA]*B[rowsB,colsB]
K2 <- Hadamard(A, B, rowsA=rowsA, rowsB=rowsB, colsA=colsA, colsB=colsB)

all.equal(K1, K2)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
A <- matrix(rnorm(30), ncol=5)
B <- matrix(rnorm(30), ncol=5)

K <- Hadamard(A, B)
c(K=pryr::address(K), A=pryr::address(A))

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
A <- Hadamard(A, B, inplace=TRUE)
c(A=pryr::address(A))  # output address remain unchanged
all.equal(K, A)        # contains the desired result 

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
A <- matrix(1:16, ncol=4)
B <- matrix(10*(1:16), ncol=4)

dimnames(A) <- list(paste("week",1:4), month.abb[1:4])
dimnames(B) <- list(c("chicken","beef","pork","fish"), LETTERS[1:4])

Hadamard(A, B, make.dimnames=TRUE)
Hadamard(A, B, make.dimnames=TRUE, colsA=c(1,1,1,1), rowsB=c(2,3,2,3))


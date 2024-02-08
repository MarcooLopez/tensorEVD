## ----initialsetup, include=FALSE----------------------------------------------
knitr::opts_chunk$set(cache=FALSE)
library(tensorEVD)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
a <- 10
( B <- matrix(1:6, ncol=2) )
Kronecker(a, B)
# In this case, should be equal to Kronecker(B, a)
Kronecker(B, a)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
D <- diag(1, 2)
Kronecker(D, B)
# Is not equal to
Kronecker(B, D)    # this a 'striped' matrix

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
a <- c(1,2,3)
b <- c(4,5)
Kronecker(a, t(b))
# Should be equal to the product a b'
tcrossprod(a, b)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Simulating matrices A and B
m = 20; n = 20
p = 40; q = 30
A <- matrix(rnorm(m*n), ncol=n)
B <- matrix(rnorm(p*q), ncol=q)

# Making the Kronecker product
K1 <- kronecker(A, B)                    
K2 <- Kronecker(A, B)                   

# Should be equal
all.equal(K1,K2)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
dm <- c(nrow(A)*nrow(B), ncol(A)*ncol(B))    # dimension of the Kronecker

# Subsetting a matrix with 30% of rows/columns
rows <- sample(seq(dm[1]), 0.3*dm[1])
cols <- sample(seq(dm[2]), 0.3*dm[2])

K1 <- Kronecker(A, B)[rows,cols]
K2 <- Kronecker(A, B, rows=rows, cols=cols)

dim(K1)   # small size
all.equal(K1, K2)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
B <- rnorm(1) # B is a scalar
K <- Kronecker(A, B)
c(K=pryr::address(K), A=pryr::address(A))

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
A <- Kronecker(A, B, inplace=TRUE)
c(A=pryr::address(A))  # output address remain unchanged
all.equal(K, A)        # contains the desired result

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
A <- matrix(1:9, ncol=3)
B <- matrix(10*(1:4), ncol=2)

dimnames(A) <- list(paste("week",1:3), month.abb[1:3])
dimnames(B) <- list(c("office","home"), LETTERS[1:2])

Kronecker(A, B, make.dimnames=TRUE)


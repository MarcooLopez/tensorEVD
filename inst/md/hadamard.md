### Hadamard product

For any two matrices
<b>A</b> = \{a<sub>ij</sub>\} and <b>B</b> = \{b<sub>ij</sub>\} of the same dimensions <i>m</i> &times; <i>n</i>, the Hadamard product between them is defined as the element-wise or entry-wise product

<p align="center">
<b>A</b> &odot; <b>B</b> = &#123;a<sub>ij</sub>b<sub>ij</sub>&#125;
</p>

This can be performed using the R's product operator ('*') 
```r
m = 20; n = 50
A = matrix(rnorm(m*n), ncol=n)
B = matrix(rnorm(m*n), ncol=n)
K = A*B
```

If <b>A</b> and <b>B</b> are not of the same dimensions (e.g., <b>A</b> is <i>m</i> &times; <i>n</i> and <b>B</b> is <i>p</i> &times; <i>q</i>), a Hadamard product can be still performed using incidence matrices of appropiate dimensions

<p align="center">
(<b>R</b><sub>1</sub> <b>A</b> <b>C'</b><sub>1</sub>) &odot; (<b>R</b><sub>2</sub> <b>B</b> <b>C'</b><sub>2</sub>)
</p>

where
<b>R</b><sub>1</sub> and <b>R</b><sub>2</sub> are incidence matrices for rows and <b>C</b><sub>1</sub> and <b>C</b><sub>2</sub> are incidence matrices for columns.

Product matrix <b>R</b><sub>1</sub> <b>A</b> <b>C'</b><sub>1</sub> can be obtained by matrix indexing using, for instance, integer vectors `rowsA` and `colsA`. Likewise, product matrix <b>R</b><sub>2</sub> <b>B</b> <b>C'</b><sub>2</sub> can be obtained using integer vectors `rowsB` and `colsB`. 
Therefore, the Hadamard product can be obtained by multiplying the indexed matrices.
For example

```r
# Simulating A and B of different dimensions
m = 20; n = 50
p = 30; q = 25
A = matrix(rnorm(m*n), ncol=n)
B = matrix(rnorm(p*q), ncol=q)

# Subsetting a matrix of 100 x 120
rowsA = sample(seq(m), 100, replace=TRUE)
rowsB = sample(seq(p), 100, replace=TRUE)
colsA = sample(seq(n), 120, replace=TRUE)
colsB = sample(seq(q), 120, replace=TRUE)

K1 = A[rowsA,colsA]*B[rowsB,colsB]
```

The *Hadamard*() function computes this Hadamard product directly from <b>A</b> and <b>B</b> without forming <b>R</b><sub>1</sub> <b>A</b> <b>C'</b><sub>1</sub> or <b>R</b><sub>2</sub> <b>B</b> <b>C'</b><sub>2</sub> matrices:

```r
K2 = Hadamard(A, B, rowsA=rowsA, rowsB=rowsB, colsA=colsA, colsB=colsB)
max(K1-K2) # should be zero
```

<- [Return](https://github.com/MarcooLopez/tensorEVD/blob/main/README.md)


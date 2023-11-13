
### Kronecker product
For any two matrices
<b>A</b> = \{a<sub>ij</sub>\} of dimensions <i>m</i> &times; <i>n</i> and
<b>B</b> = \{b<sub>ij</sub>\} of dimensions <i>p</i> &times; <i>q</i>, the direct Kronecker product between them is the matrix of dimensions <i>mp</i> &times; <i>nq</i> defined as the block matrix

<p align="center">
<b>A</b> &otimes; <b>B</b> = &#123;a<sub>ij</sub><b>B</b>&#125;
</p>

This can be computed using the *Kronecker*() function:
```r
m = 10; n = 15
p = 25; q = 20
A = matrix(rnorm(m*n), ncol=n)
B = matrix(rnorm(p*q), ncol=q)
K = Kronecker(A, B)
```

Selecting specific rows and columns from the Kronecker can be done by pre- and post- multiplication with incidence matrices 

<p align="center">
<b>R</b> (<b>A</b>&otimes;<b>B</b>) <b>C'</b>
</p>

where
<b>R</b> is an incidence matrix for rows and <b>C</b> is an incidence matrix for columns.

This sub-matrix can be obtained by indexing using, for instance, integer vectors `rows` and `cols` as:
```r
# Subsetting a matrix of 100 x 120
rows = sample(seq(m*p), 100, replace=TRUE)
cols = sample(seq(n*q), 120, replace=TRUE)
K1 = Kronecker(A, B)[rows,cols]
```

This approach computes first the Kronecker product then make the sub-setting, which can be ineffective if a relatively small number of row/columns are to be selected. The *Kronecker*() function can derive this sub-matrix directly from <b>A</b> and <b>B</b> without forming the whole Kronecker product:
```r
K2 = Kronecker(A, B, rows=rows, cols=cols)
max(K1-K2) # should be zero
```

<- [Return](https://github.com/MarcooLopez/tensorEVD/blob/main/README.md)

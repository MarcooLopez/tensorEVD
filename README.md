# Fast Factorization of High-Dimensional Tensor Product Matrices

The 'tensorEVD' R-package offers tools for calculation and factorization of high-dimensional tensor products (Hadamard and Kronecker) that are formed by smaller matrices.

*Last update: Nov 13, 2023*

## Installation and loading
Installation of 'tensorEVD' package requires a R-version greater than 3.5.0. The 'tensorEVD' package can be installed from GitHub, no CRAN version is available yet.

Installation from GitHub (developing version)
```r
install.packages('remotes',repos='https://cran.r-project.org/')  # 1. install remotes
library(remotes)                                                 # 2. load the library
install_github('MarcooLopez/tensorEVD')                          # 3. install tensorEVD from GitHub
library(tensorEVD)                                               # 4. load tensorEVD
```

## Documentation
Description of the package's main function.
```r
help(tensorEVD)  
```

## Examples
Here we present examples on the use of the functions included in the package.

* [Kronecker product](https://github.com/MarcooLopez/tensorEVD/blob/main/inst/md/kronecker.md)
* [Hadamard product](https://github.com/MarcooLopez/tensorEVD/blob/main/inst/md/hadamard.md)
* [Tensor EVD](https://github.com/MarcooLopez/tensorEVD/blob/main/inst/md/tensorEVD.md)

## Extra documentation
We provide benchmarks and an application in Genomic Prediction of the *tensorEVD*() function using data from the Genomes-to-Field (G2F) Initiative

* [Documentation](http://htmlpreview.github.io/?https://github.com/MarcooLopez/tensorEVD/blob/master/inst/doc/tensorEVD-documentation.html)

## Citation
*A manuscript is under review*

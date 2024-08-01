# Fast Factorization of High-Dimensional Tensor Product Matrices

[![CRAN status](https://www.r-pkg.org/badges/version/tensorEVD?color=green)](https://CRAN.R-project.org/package=tensorEVD)
[![CRAN checks](https://badges.cranchecks.info/worst/tensorEVD.svg)](https://cran.r-project.org/web/checks/check_results_tensorEVD.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/tensorEVD)](http://www.r-pkg.org/pkg/tensorEVD)
[![Downloads/month](https://cranlogs.r-pkg.org/badges/tensorEVD?color=blue)](http://www.r-pkg.org/pkg/tensorEVD)

The 'tensorEVD' R-package offers tools for calculation and factorization of high-dimensional tensor products (Hadamard and Kronecker) that are formed by smaller matrices.

**Funding**: [NSF PGRP-Tech Grant 2035472](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2035472&HistoricalAwards=false).

*Last update: May 22, 2024*

## Installation
Installation of 'tensorEVD' package requires an R-version &ge; 3.6.0.

From CRAN (stable version)
```r
install.packages('tensorEVD', repos='https://cran.r-project.org/')  
library(tensorEVD)                                                  
```

From GitHub (developing version)
```r
if(!'remotes' %in% rownames(installed.packages())){
  install.packages('remotes', repos='https://cran.r-project.org/')  
}
remotes::install_github('MarcooLopez/tensorEVD')                    
library(tensorEVD)                                                
```

## Documentation
Description of the package's main functions.
```r
help(package='tensorEVD', help_type='html')
```

## Examples
Here we present examples on the use of the functions included in the package.

* [Hadamard product](http://htmlpreview.github.io/?https://github.com/MarcooLopez/tensorEVD/blob/main/inst/Rmd/hadamard.html)
* [Kronecker product](http://htmlpreview.github.io/?https://github.com/MarcooLopez/tensorEVD/blob/main/inst/Rmd/kronecker.html)
* [Tensor EVD](http://htmlpreview.github.io/?https://github.com/MarcooLopez/tensorEVD/blob/main/inst/Rmd/tensorEVD.html)

## Application
We provide benchmarks and an application in Genomic Prediction of the *tensorEVD*() function using data from the Genomes-to-Field (G2F) Initiative

* Lopez-Cruz *et al.*, 2024. *G3:Genes|Genomes|Genetics* [[Manuscript](https://academic.oup.com/g3journal/article/14/3/jkae001/7511334)] [[Documentation](http://htmlpreview.github.io/?https://github.com/MarcooLopez/tensorEVD/blob/master/inst/doc/tensorEVD-documentation.html)]

## Citation
Lopez-Cruz M, Pérez-Rodríguez Paulino, and de los Campos G. **2024**. A fast algorithm to factorize high-dimensional Tensor Product matrices used in Genetic Models. *G3 Genes|Genomes|Genetics* 14(3): 1-8. doi: 10.1093/g3journal/jkae001

# idxstd4R
This package performs CPUE standardization with functionality for k-fold cross validation. This package was intially tailored for standardizing Alaska Sablefish CPUE using multiple gear types and data types, but can be genearlizable to other species and regions.

The function depends on rsample, mgcv, doMC, and a variety of other tools. To install the latest version:
```
devtools::install_github("chengmatt/idxstd4R")
```
Please see the examples folder for a general tutorial on how to implement tensor product smooths and offsets within idxstd4R. 

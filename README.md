# idxstd4R
This package performs CPUE standardization with functinality for 5-fold cross validation. This package was intially tailored for standardizing Alaska Sablefish CPUE using multiple gear types and data types, but can be genearlizable to other species and regions.

The function depends on rsample, mgcv, doMC, and a variety of other tools. To install the latest version:
```
devtools::install_github("chengmatt/idxstd4R")
```
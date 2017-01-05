# momentuHMM
An R package for Models Of animal MoveENT behavior Using multivariate Hidden Markov Models 

Get started with the vignette: [Guide to using momentuHMM](https://cran.r-project.org/web/packages/momentuHMM/vignettes/momentuHMM-guide.pdf)

## Installation instructions

### Stable release
The package is available on [CRAN](https://cran.r-project.org/web/packages/momentuHMM/index.html). To install it from CRAN,
you can use the following commands:
``` R
# install dependencies
install.packages(c("Rcpp","RcppArmadillo","sp","CircStats"))
# install momentuHMM
install.packages("momentuHMM")
```

### Install from Github
To install the latest (**unstable**) version of the package from Github:
``` R
install.packages(c("Rcpp","RcppArmadillo","sp","CircStats","devtools"))
library(devtools)
install_github("bmcclintock/momentuHMM", ref="develop", build_vignettes=TRUE)
```

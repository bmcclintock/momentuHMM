# momentuHMM [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/momentuHMM)](https://cran.r-project.org/package=momentuHMM) [![Travis-CI Build Status](https://travis-ci.org/bmcclintock/momentuHMM.svg)](https://travis-ci.org/bmcclintock/momentuHMM) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0) [![Downloads](http://cranlogs.r-pkg.org/badges/momentuHMM)](https://cran.r-project.org/package=momentuHMM)

An R package for Maximum likelihood analysis Of animal MovemENT behavior Using multivariate Hidden Markov Models 

Get started with the vignette: [Guide to using momentuHMM](https://cran.r-project.org/package=momentuHMM/vignettes/momentuHMM.pdf)

## Installation instructions

### CRAN release
The package is available on [CRAN](https://cran.r-project.org/package=momentuHMM). To install it:
``` R
install.packages("momentuHMM")
```

### Install from Github
To install the latest (stable) version of the package from Github:
``` R
library(devtools)
install_github("bmcclintock/momentuHMM")
```

To install the latest (**unstable**) version of the package from Github:
``` R
library(devtools)
install_github("bmcclintock/momentuHMM", ref="develop")
```

## References
McClintock, B.T., Michelot, T. (2018) [momentuHMM: R package for generalized hidden Markov models of animal movement](http://dx.doi.org/10.1111/2041-210X.12995). *Methods in Ecology and Evolution*. doi:10.1111/2041-210X.12995.

McClintock, B.T., King R., Thomas L., Matthiopoulos J., McConnell B.J., Morales J.M. (2012) [A general discrete-time modeling framework for animal movement using multistate random walks](http://onlinelibrary.wiley.com/doi/10.1890/11-0326.1/full). *Ecological Monographs*, 82(3), 335-349.

McClintock, B.T. (2017) [Incorporating telemetry error into hidden Markov models of animal movement using multiple imputation](https://link.springer.com/article/10.1007/s13253-017-0285-6). *Journal of Agricultural, Biological, and Environmental Statistics*, 22(3), 249-269.

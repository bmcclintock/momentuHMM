
context("trMatrix_rcpp")

test_that("The output has the right format",{
  nbCovs <- 2
  nbStates <- 3
  nbObs <- 10
  covs <- cbind(rep(1,nbObs),rnorm(nbObs),rpois(nbObs,5))
  beta <- matrix(rnorm(18),ncol=nbStates*(nbStates-1),nrow=nbCovs+1)
  trMat <- trMatrix_rcpp(nbStates,beta,covs)

  expect_equal(dim(trMat),c(nbStates,nbStates,nbObs))
})

test_that("Rows sum to 1",{
  nbCovs <- 2
  nbStates <- 3
  nbObs <- 10
  covs <- cbind(rep(1,nbObs),rnorm(nbObs),rpois(nbObs,5))
  beta <- matrix(rnorm(18),ncol=nbStates*(nbStates-1),nrow=nbCovs+1)
  trMat <- trMatrix_rcpp(nbStates,beta,covs)

  expect_equal(length(which(apply(trMat,3,rowSums)-1>1e-15)),0) # cannot expect 1 exactly
})

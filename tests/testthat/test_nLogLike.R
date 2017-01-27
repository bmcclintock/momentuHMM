
context("nLogLike")

test_that("Exceptions are thrown",{
  data <- example$data
  simPar <- example$simPar
  par0 <- example$par0

  estAngleMean <- is.null(simPar$angleMean)
  bounds <- parDef(simPar$stepDist,simPar$angleDist,simPar$nbStates,
                   estAngleMean,simPar$zeroInflation)$bounds
  parSize <- parDef(simPar$stepDist,simPar$angleDist,simPar$nbStates,
                    estAngleMean,simPar$zeroInflation)$parSize

  par <- c(par0$stepPar0,par0$anglePar0)
  wpar <- n2w(par,bounds,par0$beta0,par0$delta0,simPar$nbStates,estAngleMean)

  expect_that(nLogLike(wpar,simPar$nbStates,bounds,parSize,data,simPar$stepDist,simPar$angleDist,
                       simPar$angleMean,simPar$zeroInflation),not(throws_error()))

  # if not enough parameters provided
  expect_that(nLogLike(wpar[-1],simPar$nbStates,bounds,parSize,data,simPar$stepDist,simPar$angleDist,
                       simPar$angleMean,simPar$zeroInflation),throws_error())

  # if stepDist not in list
  expect_that(nLogLike(wpar,simPar$nbStates,bounds,parSize,data,"unif",simPar$angleDist,
                       simPar$angleMean,simPar$zeroInflation),throws_error())

  # if angleDist not in list
  expect_that(nLogLike(wpar,simPar$nbStates,bounds,parSize,data,simPar$stepDist,"norm",
                       simPar$angleMean,simPar$zeroInflation),throws_error())

  data <- data[-2] # remove data$step
  expect_that(nLogLike(wpar,simPar$nbStates,bounds,parSize,data,simPar$stepDist,simPar$angleDist,
                       simPar$angleMean,simPar$zeroInflation),throws_error())
})

test_that("angleMean=NULL, angleDist=NULL, and zeroInflation=TRUE work",{
  data <- example$data
  simPar <- example$simPar
  par0 <- example$par0

  estAngleMean <- TRUE
  bounds <- parDef(simPar$stepDist,simPar$angleDist,simPar$nbStates,
                   estAngleMean,TRUE)$bounds
  parSize <- parDef(simPar$stepDist,simPar$angleDist,simPar$nbStates,
                    estAngleMean,TRUE)$parSize

  par0$stepPar0 <- c(par0$stepPar0,rep(0.2,simPar$nbStates)) # include zero mass parameters
  par0$anglePar0 <- c(rep(0,simPar$nbStates),par0$anglePar0) # include angle mean parameters
  par <- c(par0$stepPar0,par0$anglePar0)
  wpar <- n2w(par,bounds,par0$beta0,par0$delta0,simPar$nbStates,estAngleMean)

  expect_that(nLogLike(wpar,simPar$nbStates,bounds,parSize,data,simPar$stepDist,"none",
                       NULL,TRUE),not(throws_error()))
})

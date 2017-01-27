
context("allProbs")

test_that("The output has the right format",{

  p <- allProbs(example$m,length(example$m$stateNames))

  expect_equal(nrow(p),length(example$m$data$step))
  expect_equal(ncol(p),length(example$m$stateNames))
})

test_that("It works without turning angles",{
  m<-example$m
  m$conditions$dist$angle<-NULL

  expect_error(allProbs(m,length(example$m$stateNames)),NA)
})

test_that("Zero-inflation works",{
  m<-example$m
  nbStates<-length(example$m$stateNames)
  m$conditions$zeroInflation$step<-TRUE
  m$mod$estimate<-c(m$mod$estimate[1:(2*nbStates)],boot::logit(c(0.2,0.3)),m$mod$estimate[-(1:(2*nbStates))])
  m$conditions$fullDM$step<-diag(3*nbStates)
  m$conditions$cons$step<-rep(1,3*nbStates)
  m$conditions$workcons$step<-rep(0,3*nbStates)
  m$conditions$bounds$step<-rbind(m$conditions$bounds$step,matrix(c(0,1),nbStates,2,byrow=TRUE))

  expect_error(allProbs(m,nbStates),NA)
})


context("nLogLike")

test_that("Exceptions are thrown",{
  data <- example$data
  m<-example$m
  Par <- list(step=example$par0$Par$step,angle=example$par0$Par$angle)
  nbStates <- length(m$stateNames)
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  
  wpar <- n2w(Par,m$conditions$bounds,m$mle$beta,m$mle$delta,nbStates,m$conditions$estAngleMean,m$conditions$DM,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)

  expect_error(nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,m$conditions$dist,model.matrix(m$conditions$formula,data),
                       m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,
                       m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex),NA)

})

test_that("angleDist=NULL and zeroInflation=TRUE work",{
  data <- example$data
  m<-example$m
  Par <- list(step=c(example$par0$Par$step,.05,.05))
  zeroInflation<-list(step=TRUE)
  nbStates <- length(m$stateNames)
  dist <- list(step="gamma")
  
  inputs <- checkInputs(nbStates,dist,Par,NULL,NULL,zeroInflation,NULL,NULL,NULL,NULL,m$stateNames)
  DMinputs<-getDM(data,inputs$DM,dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,inputs$cons,inputs$workcons,zeroInflation,inputs$circularAngleMean)

  wpar <- n2w(Par,inputs$p$bounds,m$mle$beta,m$mle$delta,nbStates,inputs$estAngleMean,m$conditions$DM,DMinputs$cons,DMinputs$workcons,m$conditions$Bndind)
  
  expect_error(nLogLike(wpar,nbStates,m$conditions$formula,inputs$p$bounds,inputs$p$parSize,data,dist,model.matrix(m$conditions$formula,data),
                       inputs$estAngleMean,inputs$circularAngleMean,zeroInflation,
                       m$conditions$stationary,DMinputs$cons,DMinputs$fullDM,m$conditions$DMind,DMinputs$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex),NA)
})

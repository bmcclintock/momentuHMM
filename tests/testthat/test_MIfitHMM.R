context("MIfitHMM")

test_that("Exceptions are thrown",{
  
  oldRNG<-setRNG::setRNG()
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
  
  theta<-c(6,-0.5)
  err.model <- miExample$err.model
  
  # temporally irregular data with measurement error
  obsData <- miExample$obsData
  init <- miExample$inits
  
  fixPar<-c(1,1,NA,NA)
  crwOut<-crawlWrap(obsData,ncores=1,retryFits=0,initial.state=init,err.model=err.model,theta=theta,fixPar=fixPar,attempts=100)
  
  bestData<-prepData(crwOut,covNames=names(example$m$rawCovs))
  bestFit<-fitHMM(bestData,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par0=getPar(example$m)$Par,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM)
  
  Par<-getPar(bestFit)
  
  expect_error(MIfitHMM(miData=crwOut,nSims=3,ncores=1,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par0=Par$Par,beta0=Par$beta,delta0=Par$delta,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM,covNames=names(example$m$rawCovs),method="quadrature"),NA)
  
  setRNG::setRNG(oldRNG)
  
})
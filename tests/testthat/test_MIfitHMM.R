context("MIfitHMM")

test_that("Exceptions are thrown",{
  
  oldRNG<-setRNG::setRNG()
  
  setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=2)
  
  theta<-c(6,-0.5)
  err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
  
  # temporally irregular data with measurement error
  obsData <- miExample$obsData
  
  fixPar<-c(1,1,NA,NA)
  crwOut<-crawlWrap(obsData,retryFits=1,err.model=err.model,theta=theta,fixPar=fixPar,attempts=100)
  
  bestData<-prepData(crwOut,covNames=names(example$m$rawCovs))
  bestFit<-fitHMM(bestData,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par0=getPar(example$m)$Par,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM)
  
  Par<-getPar(bestFit)
  
  expect_error(MIfitHMM(miData=crwOut,nSims=3,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par0=Par$Par,beta0=Par$beta,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM,covNames=names(example$m$rawCovs),method="quadrature"),NA)
  
  setRNG::setRNG(oldRNG)
  
})

test_that("betaCons works correctly",{
  
  oldRNG<-setRNG::setRNG()
  
  setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=3)
  
  theta<-c(6,-0.5)
  err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
  
  # temporally irregular data with measurement error
  obsData <- miExample$obsData
  
  fixPar<-c(1,1,NA,NA)
  crwOut<-crawlWrap(obsData,retryFits=1,err.model=err.model,theta=theta,fixPar=fixPar,attempts=100)
  
  bestData<-prepData(crwOut,covNames=names(example$m$rawCovs))
  bestFit<-fitHMM(bestData,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par0=getPar(example$m)$Par,formula=example$m$conditions$formula,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM,betaCons=matrix(c(1,2,3,4,2,3),3,2))
  
  Par<-getPar(bestFit)
  
  expect_error(MIfitHMM(miData=crwOut,nSims=3,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par0=Par$Par,beta0=Par$beta,formula=example$m$conditions$formula,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM,betaCons=matrix(c(1,2,3,4,2,3),3,2),covNames=names(example$m$rawCovs),method="quadrature"),NA)
  
  setRNG::setRNG(oldRNG)
  
})

context("doFuture")

test_that("Parallel processing works",{
  obsData <- miExample$obsData
  
  err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
  
  ncores <- 2
  
  set.seed(6,kind="Mersenne-Twister",normal.kind = "Inversion")
  expect_error(crwOut <- crawlWrap(obsData=obsData,ncores=ncores,retryFits=10,
                      theta=c(4,0),fixPar=c(1,1,NA,NA),
                      err.model=err.model,attempts=100),NA)
  
  bPar <- miExample$bPar
  expect_error(HMMfits <- MIfitHMM(crwOut,nSims=4,ncores=ncores,
                      nbStates=2,dist=list(step="gamma",angle="vm"),
                      Par0=bPar$Par,beta0=bPar$beta,
                      formula=~cov1+cos(cov2),
                      estAngleMean=list(angle=TRUE),
                      covNames=c("cov1","cov2")),NA)
  
  expect_error(plotPR(HMMfits,ncores=ncores),NA)
  
  expect_error(timeInStates(HMMfits$HMMfits,ncores=ncores),NA)
  
  expect_error(simData(model=HMMfits,nbAnimals=2,obsPerAnimal=30,ncores=ncores),NA)
})
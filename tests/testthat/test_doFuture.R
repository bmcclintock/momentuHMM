
context("doFuture")

test_that("Parallel processing works",{
  
  ncores <- 2
  
  nbStates <- 2
  stepDist <- "gamma" # step distribution
  angleDist <- "vm" # turning angle distribution
  nbStates <- 2
  mu0 <- c(20,70)
  sigma0 <- c(10,30)
  kappa0 <- c(1,1)
  stepPar <- c(mu0,sigma0) # no zero-inflation, so no zero-mass included
  anglePar <- kappa0 # not estimating angle mean, so not included
  formula <- ~cov1+cos(cov2)
  
  expect_error(m <- fitHMM(data=example$m$data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
              Par0=list(step=stepPar,angle=anglePar),formula=formula,retryFits=4,ncores=ncores,nlmPar=list(hessian=FALSE)),NA)
  
  obsData <- miExample$obsData
  
  err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
  
  set.seed(2,kind="Mersenne-Twister",normal.kind = "Inversion")
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
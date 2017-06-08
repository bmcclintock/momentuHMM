context("crawlWrap")

test_that("Exceptions are thrown",{

  oldRNG<-setRNG::setRNG()
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
  
  theta<-c(6,-0.5)
  err.model <- miExample$err.model
  
  # temporally irregular data
  obsData <- simData(model=example$m,lambda=2)
  init <- list(a = c(obsData$x[1],0,
                     obsData$y[1],0),
               P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                          5000 ^ 2, 10 * 3600 ^ 2)))
  fixPar<-c(NA,NA)
  expect_error(crawlWrap(obsData,ncores=1,initial.state=init,theta=theta,fixPar=fixPar,attempts=100,retryFits=0),NA)
  
  # measurement error
  obsData <- simData(model=example$m,errorEllipse=list(M=50,m=50,r=0))
  init <- list(a = c(obsData$x[1],0,
                     obsData$y[1],0),
               P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                          5000 ^ 2, 10 * 3600 ^ 2)))
  fixPar<-c(1,1,NA,NA)
  expect_error(crawlWrap(obsData,ncores=1,initial.state=init,err.model=err.model,theta=theta,fixPar=fixPar,attempts=100,retryFits=0),NA)
  
  # temporally irregular data with measurement error
  obsData <- miExample$obsData
  init <- miExample$inits
  
  fixPar<-c(1,1,NA,NA)
  expect_error(crawlWrap(obsData,ncores=1,initial.state=init,err.model=err.model,theta=theta,fixPar=fixPar,attempts=100,retryFits=0),NA)
  
  # temporally regular data without measurement error data
  obsData <- simData(model=example$m)
  fixPar<-c(NA,NA)
  expect_error(crawlWrap(obsData,ncores=1,initial.state=init,theta=theta,fixPar=fixPar))
  
  setRNG::setRNG(oldRNG)
  
})
context("crawlWrap")

test_that("Exceptions are thrown",{

  oldRNG<-setRNG::setRNG()
  
  setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=1)
  
  theta<-c(6,-0.5)
  err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
  
  # temporally irregular data
  obsData <- simData(model=example$m,lambda=2)
  fixPar<-c(NA,NA)
  expect_error(crawlWrap(obsData,theta=theta,fixPar=fixPar,attempts=100,retryFits=0),NA)
  
  # measurement error
  obsData <- simData(model=example$m,errorEllipse=list(M=50,m=50,r=0))

  fixPar<-c(1,1,NA,NA)

  expect_error(crwOut<-crawlWrap(obsData,ncores=1,err.model=err.model,theta=theta,fixPar=fixPar,attempts=100,retryFits=0),NA)
  expect_error(plot(crwOut,ask=FALSE),NA)
  expect_error(plot(crwOut,ask=FALSE,crawlPlot=TRUE),NA)
  
  # temporally irregular data with measurement error
  obsData <- miExample$obsData
  
  fixPar<-c(1,1,NA,NA)

  expect_error(crwOut<-crawlWrap(obsData,ncores=1,err.model=err.model,theta=theta,fixPar=fixPar,attempts=100,retryFits=0),NA)
  expect_error(plot(crwOut,ask=FALSE),NA)
  expect_error(plot(crwOut,ask=FALSE,crawlPlot=TRUE),NA)
  
  # temporally regular data without measurement error data
  obsData <- simData(model=example$m)
  fixPar<-c(NA,NA)
  expect_error(crawlWrap(obsData,theta=theta,fixPar=fixPar))
  
  setRNG::setRNG(oldRNG)
  
})
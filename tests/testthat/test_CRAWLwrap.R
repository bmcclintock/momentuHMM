context("CRAWLwrap")

test_that("Exceptions are thrown",{

  theta<-c(6,-0.5)
  err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
  
  # temporally irregular data
  obsData <- simData(model=example$m,lambda=1)
  init <- list(a = c(obsData$x[1],0,
                     obsData$y[1],0),
               P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                          5000 ^ 2, 10 * 3600 ^ 2)))
  fixPar<-c(NA,NA)
  expect_error(CRAWLwrap(obsData,ncores=1,initial.state=init,theta=theta,fixPar=fixPar),NA)
  
  # measurement error
  obsData <- simData(model=example$m,errorEllipse=list(M=150,m=150,r=180))
  init <- list(a = c(obsData$x[1],0,
                     obsData$y[1],0),
               P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                          5000 ^ 2, 10 * 3600 ^ 2)))
  fixPar<-c(1,1,NA,NA)
  expect_error(CRAWLwrap(obsData,ncores=1,initial.state=init,err.model=err.model,theta=theta,fixPar=fixPar),NA)
  
  # temporally irregular data with measurement error
  obsData <- simData(model=example$m,lambda=0.5,errorEllipse=list(M=150,m=150,r=180))
  init <- list(a = c(obsData$x[1],0,
                     obsData$y[1],0),
               P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                          5000 ^ 2, 10 * 3600 ^ 2)))
  fixPar<-c(1,1,NA,NA)
  expect_error(CRAWLwrap(obsData,ncores=1,initial.state=init,err.model=err.model,theta=theta,fixPar=fixPar),NA)
  
  # temprally regular data without mwasurement error data
  fixPar<-c(NA,NA)
  expect_error(CRAWLwrap(simData(model=example$m),ncores=1,initial.state=init,theta=theta,fixPar=fixPar))
  
})
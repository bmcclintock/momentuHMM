context("MIfitHMM")

test_that("Exceptions are thrown",{
  
  theta<-c(6,-0.5)
  err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
  
  # temporally irregular data with measurement error
  obsData <- simData(model=example$m,lambda=2,errorEllipse=list(M=10,m=10,r=180))
  init <- list(a = c(obsData$mux[1],0,
                     obsData$muy[1],0),
               P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                          5000 ^ 2, 10 * 3600 ^ 2)))
  fixPar<-c(1,1,NA,NA)
  crwOut<-CRAWLwrap(obsData,ncores=1,retryFits=100,initial.state=init,err.model=err.model,theta=theta,fixPar=fixPar,attempts=20)
  
  
  bestFit<-MIfitHMM(1,ncores=1,miData=crwOut,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par=getPar(example$m)$Par,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM,covNames=names(example$m$rawCovs))
  
  Par<-getPar(bestFit)
  
  expect_error(MIfitHMM(2,ncores=1,miData=crwOut,nbStates=length(example$m$stateNames),dist=example$m$conditions$dist,Par=Par$Par,beta0=Par$beta,delta0=Par$delta,estAngleMean=example$m$conditions$estAngleMean,DM=example$m$conditions$DM,covNames=names(example$m$rawCovs)),NA)
  
})

context("getDM")

test_that("Exception is thrown",{
  
  m<-example$m
  p<-parDef(m$conditions$dist,length(m$stateNames),m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds)
  expect_error(getDM(m$data,m$conditions$DM,m$conditions$dist,length(m$stateNames),p$parNames,p$bounds,getPar(m)$Par,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean,ParChecks=TRUE),
              NA)
})

test_that("List DM and matrix DM are the same",{

nbStates <- 2
stepDist <- "gamma" 
angleDist <- "vm" 

data <- example$m$data

mu0 <- c(20,70)
sigma0 <- c(10,30)
kappa0 <- c(1,1)
stepPar <- c(mu0,sigma0) 
anglePar <- kappa0 
formula <- ~cov1+cos(cov2)

listDM<-list(step=list(mean=formula,sd=~1))
matrixDM<-list(step=matrix(c(1,"cov1","cos(cov2)",0,0,0,0,0,
                             0,0,0,1,"cov1","cos(cov2)",0,0,
                             0,0,0,0,0,0,1,0,
                             0,0,0,0,0,0,0,1),nrow=2*nbStates,byrow=TRUE,dimnames = list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates)),c(paste0("mean_",rep(1:nbStates,each=3),":",c("(Intercept)","cov1","cos(cov2)")),paste0("sd_",1:nbStates,":(Intercept)")))))

listPar0 <- getParDM(data,nbStates,list(step=stepDist,angle=angleDist),
                 Par=list(step=stepPar,angle=anglePar),
                 DM=listDM)
listmod<-fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
                Par0=listPar0,
                formula=formula,
                DM=listDM)

matrixPar0 <- getParDM(data,nbStates,list(step=stepDist,angle=angleDist),
                     Par=list(step=stepPar,angle=anglePar),
                     DM=matrixDM)
matrixmod<-fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
                  Par0=matrixPar0,
                  formula=formula,
                  DM=matrixDM)

  expect_equal(all.equal(listmod$conditions$fullDM,matrixmod$conditions$fullDM),TRUE)
  expect_equal(listmod$mod$minimum,matrixmod$mod$minimum)
  expect_equal(all.equal(listmod$mod$estimate,matrixmod$mod$estimate),TRUE)
})

test_that("NULL, list, and matrix DM are the same",{
  
  nbStates <- 2
  stepDist <- "gamma" 
  angleDist <- "vm" 
  
  data <- example$m$data
  
  mu0 <- c(20,70)
  sigma0 <- c(10,30)
  kappa0 <- c(1,1)
  stepPar <- c(mu0,sigma0) 
  anglePar <- kappa0 
  formula <- ~cov1+cos(cov2)
  
  listDM<-list(step=list(mean=~1,sd=~1))
  matrixDM<-list(step=matrix(diag(2*nbStates),nrow=2*nbStates,byrow=TRUE,dimnames = list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates)),c(paste0("mean_",1:nbStates,":(Intercept)"),paste0("sd_",1:nbStates,":(Intercept)")))))
  
  nullmod<-fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
                  Par0=list(step=stepPar,angle=anglePar),beta0=example$m$mle$beta,delta0=example$m$mle$delta[1,],
                  formula=formula,fit=FALSE)
  
  listPar0 <- getParDM(data,nbStates,list(step=stepDist,angle=angleDist),
                       Par=list(step=stepPar,angle=anglePar),
                       DM=listDM)
  listmod<-fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
                  Par0=listPar0,beta0=example$m$mle$beta,delta0=example$m$mle$delta[1,],
                  formula=formula,
                  DM=listDM,fit=FALSE)
  
  matrixPar0 <- getParDM(data,nbStates,list(step=stepDist,angle=angleDist),
                         Par=list(step=stepPar,angle=anglePar),
                         DM=matrixDM)
  matrixmod<-fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
                    Par0=matrixPar0,beta0=example$m$mle$beta,delta0=example$m$mle$delta[1,],
                    formula=formula,
                    DM=matrixDM,fit=FALSE)
  
  expect_equal(all.equal(listmod$conditions$fullDM,matrixmod$conditions$fullDM),all.equal(listmod$conditions$fullDM,nullmod$conditions$fullDM))
})

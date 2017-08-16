
context("nLogLike")

test_that("Exceptions are thrown",{
  data <- example$m$data
  m<-example$m
  Par <- list(step=example$par0$Par$step,angle=example$par0$Par$angle)
  nbStates <- length(m$stateNames)
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  
  distnames<-names(m$conditions$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
  }
  
  wpar <- n2w(Par,m$conditions$bounds,m$mle$beta,m$mle$delta,nbStates,m$conditions$estAngleMean,m$conditions$DM,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)

  expect_error(nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,m$conditions$dist,model.matrix(m$conditions$formula,data),
                       m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,
                       m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind),NA)

})

test_that("logAlpha, logBeta, and nLogLike are consistent",{
  data <- example$m$data
  m<-example$m
  Par <- getPar(m)$Par
  nbStates <- length(m$stateNames)
  nbAnimals<-length(unique(data$ID))
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  
  distnames<-names(m$conditions$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
  }
  
  wpar <- n2w(Par,m$conditions$bounds,m$mle$beta,m$mle$delta,nbStates,m$conditions$estAngleMean,m$conditions$DM,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)

  # all data
  nll<-nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,m$conditions$dist,model.matrix(m$conditions$formula,data),
                       m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,
                       m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind)
  la<-logAlpha(m)
  lb<-logBeta(m)
  ll<-0
  for(i in 1:nbAnimals){
    aInd<-max(which(data$ID==i))
    c <- max(la[aInd,]+lb[aInd,]) # cancels below ; prevents numerical errors
    ll <- ll + c + log(sum(exp(la[aInd,]+lb[aInd,]-c)))
  }
  expect_equal(nll,-ll)
  
  # random time step from each individual
  for(i in 1:nbAnimals){
    data<-example$m$data
    aInd<-which(data$ID==i)
    
    data<-data[aInd,]
    nll<-nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,m$conditions$dist,model.matrix(m$conditions$formula,data),
                       m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,
                       m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind)
  
    samp<-sample(aInd,1)
    c <- max(la[samp,]+lb[samp,]) # cancels below ; prevents numerical errors
    ll <- c + log(sum(exp(la[samp,]+lb[samp,]-c)))
    expect_equal(nll,-ll)
  }
  
})

test_that("angleDist=NULL and zeroInflation=TRUE work",{
  data <- example$m$data
  m<-example$m
  Par <- list(step=c(example$par0$Par$step,.05,.05))
  zeroInflation<-list(step=TRUE)
  oneInflation<-list(step=FALSE)
  nbStates <- length(m$stateNames)
  dist <- list(step="gamma")
  
  inputs <- checkInputs(nbStates,dist,Par,NULL,NULL,zeroInflation,oneInflation,NULL,NULL,NULL,NULL,m$stateNames)
  DMinputs<-getDM(data,inputs$DM,dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,inputs$cons,inputs$workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
  
  distnames<-names(m$conditions$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
  }

  wpar <- n2w(Par,inputs$p$bounds,m$mle$beta,m$mle$delta,nbStates,inputs$estAngleMean,m$conditions$DM,DMinputs$cons,DMinputs$workcons,m$conditions$Bndind)
  
  expect_error(nLogLike(wpar,nbStates,m$conditions$formula,inputs$p$bounds,inputs$p$parSize,data,dist,model.matrix(m$conditions$formula,data),
                       inputs$estAngleMean,inputs$circularAngleMean,zeroInflation,oneInflation,
                       m$conditions$stationary,DMinputs$cons,DMinputs$fullDM,m$conditions$DMind,DMinputs$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind),NA)
})


context("nLogLike")

test_that("Exceptions are thrown",{
  data <- example$m$data
  m<-example$m
  Par <- list(step=example$par0$Par$step,angle=example$par0$Par$angle)
  nbStates <- length(m$stateNames)
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  
  distnames<-names(inputs$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) {
      meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  
  wpar <- n2w(Par,m$conditions$bounds,m$mle$beta,log(m$mle$delta[-1]/m$mle$delta[1]),nbStates,m$conditions$estAngleMean,m$conditions$DM,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)

  expect_error(nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,inputs$dist,model.matrix(m$conditions$formula,data),
                       m$conditions$estAngleMean,m$conditions$circularAngleMean,inputs$consensus,m$conditions$zeroInflation,m$conditions$oneInflation,
                       m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind,m$covsDelta,m$conditions$workBounds),NA)

})

test_that("logAlpha, logBeta, and nLogLike are consistent",{
  data <- example$m$data
  m<-example$m
  Par <- getPar(m)$Par
  nbStates <- length(m$stateNames)
  nbAnimals<-length(unique(data$ID))
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  
  distnames<-names(inputs$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) {
      meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  
  nbCovsDelta <- ncol(m$covsDelta)-1
  foo <- length(m$mod$estimate)-(nbCovsDelta+1)*(nbStates-1)+1
  delta<- m$mod$estimate[foo:length(m$mod$estimate)]
  
  wpar <- n2w(Par,m$conditions$bounds,m$mle$beta,delta,nbStates,m$conditions$estAngleMean,m$conditions$DM,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)

  # all data
  nll<-nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,inputs$dist,model.matrix(m$conditions$formula,data),
                       m$conditions$estAngleMean,m$conditions$circularAngleMean,inputs$consensus,m$conditions$zeroInflation,m$conditions$oneInflation,
                       m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind,m$covsDelta,m$conditions$workBounds)
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
    nll<-nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,inputs$dist,model.matrix(m$conditions$formula,data),
                       m$conditions$estAngleMean,m$conditions$circularAngleMean,inputs$consensus,m$conditions$zeroInflation,m$conditions$oneInflation,
                       m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind,m$covsDelta,m$conditions$workBounds)
  
    samp<-sample(aInd,1)
    c <- max(la[samp,]+lb[samp,]) # cancels below ; prevents numerical errors
    ll <- c + log(sum(exp(la[samp,]+lb[samp,]-c)))
    expect_equal(nll,-ll)
  }
  
})

test_that("logAlpha, logBeta, and nLogLike are consistent with zero and one inflation",{
  
  oldRNG<-setRNG::setRNG()
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
  
  data <- simData(nbStates=2,dist=list(dive="beta"),Par=list(dive=c(1,1,10,1,0.2,0.1,0.1,0.2)),zeroInflation = list(dive=TRUE),oneInflation = list(dive=TRUE))
  m<-fitHMM(data,nbStates=2,dist=list(dive="beta"),Par0=list(dive=c(1,1,10,1,0.2,0.1,0.1,0.2)))
  Par <- getPar(m)$Par
  nbStates <- length(m$stateNames)
  nbAnimals<-length(unique(data$ID))
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  
  distnames<-names(inputs$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) {
      meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  
  nbCovsDelta <- ncol(m$covsDelta)-1
  foo <- length(m$mod$estimate)-(nbCovsDelta+1)*(nbStates-1)+1
  delta<- m$mod$estimate[foo:length(m$mod$estimate)]
  
  wpar <- n2w(Par,m$conditions$bounds,m$mle$beta,delta,nbStates,m$conditions$estAngleMean,m$conditions$DM,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)
  
  # all data
  nll<-nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,inputs$p$parSize,data,inputs$dist,model.matrix(m$conditions$formula,data),
                m$conditions$estAngleMean,m$conditions$circularAngleMean,inputs$consensus,m$conditions$zeroInflation,m$conditions$oneInflation,
                m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind,m$covsDelta,m$conditions$workBounds)
  la<-logAlpha(m)
  lb<-logBeta(m)
  ll<-0
  for(i in 1:nbAnimals){
    aInd<-max(which(data$ID==i))
    c <- max(la[aInd,]+lb[aInd,]) # cancels below ; prevents numerical errors
    ll <- ll + c + log(sum(exp(la[aInd,]+lb[aInd,]-c)))
  }
  expect_equal(nll,-ll)
  
  setRNG::setRNG(oldRNG)
  
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
  DMinputs<-getDM(data,inputs$DM,inputs$dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,inputs$cons,inputs$workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
  
  distnames<-names(inputs$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) {
      meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }

  wpar <- n2w(Par,inputs$p$bounds,m$mle$beta,m$mle$delta,nbStates,inputs$estAngleMean,m$conditions$DM,DMinputs$cons,DMinputs$workcons,m$conditions$Bndind)
  
  expect_error(nLogLike(wpar,nbStates,m$conditions$formula,inputs$p$bounds,inputs$p$parSize,data,inputs$dist,model.matrix(m$conditions$formula,data),
                       inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,zeroInflation,oneInflation,
                       m$conditions$stationary,DMinputs$cons,DMinputs$fullDM,m$conditions$DMind,DMinputs$workcons,m$conditions$Bndind,m$knownStates,m$conditions$fixPar,m$conditions$wparIndex,nc,meanind,m$covsDelta,m$conditions$workBounds),NA)
})

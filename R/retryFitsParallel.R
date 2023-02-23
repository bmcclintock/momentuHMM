retryFitsParallel <- function(data,nbStates, dist, Par0, beta0, delta0,
                              estAngleMean, circularAngleMean, formula, formulaDelta, stationary, mixtures, formulaPi,
                              nlmPar, fit, DM,
                              userBounds, workBounds, betaCons, betaRef, deltaCons, mvnCoords, stateNames, knownStates, fixPar, retryFits, retrySD, ncores, optMethod, control, prior, modelName, ...){
  
  for(pkg in c("doFuture","future")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package \"",pkg,"\" needed for parallel processing to work. Please install it.",
           call. = FALSE)
    }
  }
  oldDoPar <- doFuture::registerDoFuture()
  on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
  future::plan(future::multisession, workers = ncores)
  # hack so that foreach %dorng% can find internal momentuHMM variables without using ::: (forbidden by CRAN)
  progBar <- progBar
  pkgs <- c("momentuHMM")
  quietCrawl <- quietCrawl
  
  if(optMethod=="nlm"){
    hess <- nlmPar$hessian
    tmpPar <- nlmPar
    tmpPar$hessian <- FALSE
    nlmPar <- tmpPar
    control <- NULL
  } else {
    hess <- control$hessian
    tmpPar <- control
    tmpPar$hessian <- FALSE
    control <- tmpPar
    nlmPar <- NULL
  }
  curmod <- fitHMM(data,nbStates, dist, Par0, beta0, delta0,
                   estAngleMean, circularAngleMean, formula, formulaDelta, stationary, mixtures, formulaPi,
                   nlmPar=nlmPar, fit, DM,
                   userBounds, workBounds, betaCons, betaRef, deltaCons, mvnCoords, stateNames, knownStates, fixPar, retryFits=0, retrySD, ncores=1, optMethod, control=control, prior, modelName, ...)
  
  cat("\nInitial log-likelihood value: ",-curmod$mod$minimum,"\n",sep="")
  
  message("\nAttempting to improve fit in parallel using ",retryFits," random perturbations. Press 'esc' to force exit from 'fitHMM'")
  
  curPar <- getPar0(curmod)
  
  progBar(0,retryFits)
  
  dontcheck <- ifelse(isTRUE(list(...)$CT),TRUE,FALSE)
  
  
  j <- NULL
  withCallingHandlers(fits <- foreach(j = 1:retryFits, .export=c("fitHMM","fitHMM.momentuHMMData","fitHMM.momentuHierHMMData"), .errorhandling="pass", .packages=pkgs, .inorder = FALSE) %dorng% {
                                 
                                 progBar(j,retryFits)
                                 
                                 quietCrawl(jFit<-suppressWarnings(suppressMessages(fitHMM(data=data,nbStates=nbStates,dist=dist,
                                                                          Par0=curPar$Par,beta0=curPar$beta,delta0=curPar$delta,
                                                                          estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                                                                          formula=formula,formulaDelta=formulaDelta,stationary=stationary,mixtures=mixtures,formulaPi=formulaPi,
                                                                          nlmPar=nlmPar,fit=fit,
                                                                          DM=DM,userBounds=userBounds,workBounds=workBounds,betaCons=betaCons,betaRef=betaRef,deltaCons=deltaCons,
                                                                          mvnCoords=mvnCoords,stateNames=stateNames,knownStates=knownStates,fixPar=fixPar,retryFits=1,retrySD=retrySD,ncores=1,optMethod=optMethod,control=control,prior=prior,modelName=modelName, ..., dontcheck=dontcheck))))
                                 jFit
                               } 
                             ,warning=muffleCTwarning)
  future::plan(future::sequential)
  if(!all(unlist(lapply(fits,function(x) inherits(x,"error"))))){
    bestInd <- which.min(unlist(lapply(fits,function(x) ifelse(!inherits(x,"error"),x$mod$minimum,Inf))))
    curmod <- fits[[bestInd]]
    if(!isFALSE(hess)){
      cat("Calculating hessian for best model fit...\n\n")
      curPar <- getPar(curmod)
      if(optMethod=="nlm"){
        nlmPar$hessian <- TRUE
        nlmPar$print.level <- 0
      } else {
        control$hessian <- TRUE
        control$trace <- 0
      }
      curmod <- suppressMessages(fitHMM(data,nbStates, dist, curPar$Par, curPar$beta, curPar$delta,
                       estAngleMean, circularAngleMean, formula, formulaDelta, stationary, mixtures, formulaPi,
                       nlmPar=nlmPar, fit, DM,
                       userBounds, workBounds, betaCons, betaRef, deltaCons, mvnCoords, stateNames, knownStates, fixPar, retryFits=0, retrySD, ncores=1, optMethod, control=control, prior, modelName, ...))
    }
  } else warning("all retryFits failed")
  cat("Final log-likelihood value: ",-curmod$mod$minimum,"\n",sep="")
  return(curmod)
}
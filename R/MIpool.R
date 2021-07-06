#'
#' Calculate pooled parameter estimates and states across multiple imputations
#' 
#' @param im List comprised of \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} objects
#' @param alpha Significance level for calculating confidence intervals of pooled estimates (including location error ellipses). Default: 0.95.
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing).
#' @param covs Data frame consisting of a single row indicating the covariate values to be used in the calculation of pooled natural parameters. 
#' For any covariates that are not specified using \code{covs}, the means of the covariate(s) across the imputations are used 
#' (unless the covariate is a factor, in which case the first factor in the data is used). By default, no covariates are specified.
#' @param na.rm Logical indicating whether or not to exclude model fits with \code{NA} parameter estimates or standard errors from pooling. Default: FALSE.
#' 
#' @return A \code{\link{miSum}} object, i.e., a list comprised of model and pooled parameter summaries, including \code{data} (averaged across imputations), \code{conditions}, \code{Par}, and \code{MIcombine} 
#' (as returned by \code{\link[mitools]{MIcombine}} for working parameters).
#' 
#' \code{miSum$Par} is a list comprised of:
#' \item{beta}{Pooled estimates for the working parameters}
#' \item{real}{Estimates for the natural parameters based on pooled working parameters and covariate means (or \code{covs}) across imputations (if applicable)}
#' \item{timeInStates}{The proportion of time steps assigned to each state}
#' \item{states}{The most freqent state assignment for each time step based on the \code{\link{viterbi}} algorithm for each model fit}
#' \item{stateProbs}{Pooled state probability estimates for each time step}
#' \item{mixtureProbs}{Pooled mixture probabilities for each individual (only applies if \code{mixtures>1})}
#' \item{hierStateProbs}{Pooled state probability estimates for each time step at each level of the hierarchy (only applies if \code{im} is comprised of \code{\link{momentuHierHMM}} objects)}
#' 
#' @details
#' Pooled estimates, standard errors, and confidence intervals are calculated using standard multiple imputation formulas. Working scale parameters are pooled
#' using \code{\link[mitools]{MIcombine}} and t-distributed confidence intervals. Natural scale parameters and normally-distributed confidence intervals are calculated by transforming the pooled working scale parameters 
#' and, if applicable, are based on covariate means across all imputations (and/or values specified in \code{covs}).
#' 
#' The calculation of pooled error ellipses uses \code{\link[car]{dataEllipse}} from the \code{car} package. The suggested package \code{car} is not automatically imported by \code{momentuHMM} and must be installed in order to calculate error ellipses. A warning will be triggered if the \code{car} package is required but not installed.
#' 
#' Note that pooled estimates for \code{timeInStates} and \code{stateProbs} do not include within-model uncertainty and are based entirely on across-model variability.
#' 
#' @examples
#' \dontshow{
#' set.seed(3,kind="Mersenne-Twister",normal.kind="Inversion")
#' }
#' \dontrun{
#' # Extract data and crawl inputs from miExample
#' obsData <- miExample$obsData
#' 
#' # error ellipse model
#' err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
#' 
#' # Fit crawl to obsData
#' crwOut <- crawlWrap(obsData,theta=c(4,0),fixPar=c(1,1,NA,NA),
#'                     err.model=err.model)
#'                     
#' # Fit four imputations
#' bPar <- miExample$bPar
#' HMMfits <- MIfitHMM(crwOut,nSims=4,poolEstimates=FALSE,
#'                    nbStates=2,dist=list(step="gamma",angle="vm"),
#'                    Par0=bPar$Par,beta0=bPar$beta,
#'                    formula=~cov1+cos(cov2),
#'                    estAngleMean=list(angle=TRUE),
#'                    covNames=c("cov1","cov2"))
#'                    
#' # Pool estimates
#' miSum <- MIpool(HMMfits)
#' print(miSum)
#' }
#' @export
#' @importFrom stats median var qt
#' @importFrom CircStats circ.mean
# #' @importFrom car dataEllipse
# #' @importFrom mitools MIcombine
MIpool<-function(im, alpha=0.95, ncores=1, covs=NULL, na.rm=FALSE){
  
  goodIndex <- 1:length(im)
  simind <- which((unlist(lapply(im,is.momentuHMM))))
  nsims <- length(simind)
  if(nsims<1) stop("'HMMfits' must be a list comprised of momentuHMM objects")
  
  if(ncores>1){
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
  } else { 
    doParallel::registerDoParallel(cores=ncores)
    pkgs <- NULL
  }
  
  if (!requireNamespace("mitools", quietly = TRUE)) {
    stop("Package \"mitools\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  checkmove <- which(!(unlist(lapply(im,is.momentuHMM))))
  if(length(checkmove)) {
    im[checkmove]<-NULL
    warning("The following imputations are not momentuHMM objects and will be ignored: ",paste(checkmove,collapse=", "))
    goodIndex <- goodIndex[-checkmove]
  }
  
  # check modelName
  checkNames <- lapply(im,function(x) x[match("modelName",names(x))])
  if(any(!unlist(lapply(checkNames,function(x) isTRUE(all.equal(x,checkNames[[1]],use.names=FALSE)))))) stop("'modelName' must be identical for each fitted model")
  
  checksims <- lapply(im,function(x) x[match("conditions",names(x))])
  ident <- !unlist(lapply(checksims,function(x) isTRUE(all.equal(x,checksims[[1]]))))
  if(any(ident)){
    # check that only differences are in the design matrix covariate values
    checksims2 <- lapply(checksims, function(x) x$conditions[-match(c("fullDM","hierBeta","hierDelta"),names(x$conditions),nomatch=0)])
    ident2 <- !unlist(lapply(checksims2,function(x) isTRUE(all.equal(x,checksims2[[1]]))))
    if(any(ident2)) stop("Model conditions for each imputation must be identical. Imputations that do not match the first: ",paste(which(ident),collapse=", "))
  }
  
  if(any(unlist(lapply(im,function(x) is.null(x$mod$hessian))))) stop("Estimates cannot be pooled unless Hessian is calculated. Hessian is missing for imputations ",paste0(which(unlist(lapply(im,function(x) is.null(x$mod$hessian)))),collapse=", "))
  
  tmpDet <- which(unlist(lapply(im,function(x) det(x$mod$hessian)))==0)
  if(length(tmpDet)){
    warning("Hessian is singular for HMM fit(s): ",paste0(goodIndex[tmpDet],collapse=", "))
  }
  
  tmpVar <- which(unlist(lapply(im,function(x) inherits(x$mod$Sigma,"error"))))
  if(length(tmpVar)){
    warning("ginv of the hessian failed for HMM fit(s): ",paste0(goodIndex[tmpVar],collapse=", "))
    im[tmpVar] <- NULL
    nsims <- length(im)
    if(nsims<2) stop("Pooling requires at least 2 valid HMM fits")
    goodIndex <- goodIndex[-tmpVar]
  }
  
  im <- lapply(im,delta_bc)
  m <- im[[1]]
  
  wBounds <- cbind(unlist(lapply(m$conditions$workBounds,function(x) x[,1])),unlist(lapply(m$conditions$workBounds,function(x) x[,2])))
  
  # check for finite coefficients and standard errors
  betaVar <- lapply(im,function(x) get_gradwb(x$mod$estimate,wBounds)%*%x$mod$Sigma%*%t(get_gradwb(x$mod$estimate,wBounds)))
  betaCoeff <- lapply(im,function(x) w2wn(x$mod$estimate,wBounds))
  tmpVar1 <- which(unlist(lapply(betaCoeff,function(x) any(!is.finite(x)))))
  if(length(tmpVar1)){
    warning("working parameter estimates are not finite for HMM fits ",paste0(goodIndex[tmpVar1],collapse=", "),ifelse(na.rm," and will not be included in pooling",""))
    if(na.rm){
      im[tmpVar1] <- NULL
      nsims <- length(im)
      if(nsims<2) stop("Pooling requires at least 2 valid HMM fits")
      m <- im[[1]]
      betaVar <- lapply(im,function(x) get_gradwb(x$mod$estimate,wBounds)%*%x$mod$Sigma%*%t(get_gradwb(x$mod$estimate,wBounds)))
      betaCoeff <- lapply(im,function(x) w2wn(x$mod$estimate,wBounds))
      goodIndex <- goodIndex[-tmpVar1]
    } else {
      betaCoeff[tmpVar1] <- lapply(betaCoeff[tmpVar1],function(x) {x[which(!is.finite(x))]<-NA; return(x)})
    }
  }
  tmpVar2 <- unique(c(which(unlist(lapply(betaVar,function(x) any(!is.finite(x))))),which(unlist(lapply(betaVar,function(x) any(!is.finite(sqrt(diag(x)))))))))
  if(length(tmpVar2)){
    warning("working parameter standard errors are not finite for HMM fits ",paste0(goodIndex[tmpVar2],collapse=", "),ifelse(na.rm," and will not be included in pooling",""))
    if(na.rm){
      im[tmpVar2] <- NULL
      nsims <- length(im)
      if(nsims<2) stop("Pooling requires at least 2 valid HMM fits")
      m <- im[[1]]
      betaCoeff <- lapply(im,function(x) w2wn(x$mod$estimate,wBounds))
      betaVar <- lapply(im,function(x) get_gradwb(x$mod$estimate,wBounds)%*%x$mod$Sigma%*%t(get_gradwb(x$mod$estimate,wBounds)))
      goodIndex <- goodIndex[-tmpVar2]
    } else {
      betaVar[tmpVar2] <- suppressWarnings(lapply(betaVar[tmpVar2],function(x) {x[which(!is.finite(x))]<-NA; diag(x)[which(!is.finite(sqrt(diag(x))))]<-NA; return(x)}))
    }
  }
  
  data <- m$data
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))
  dist <- m$conditions$dist
  distnames <- names(dist)
  estAngleMean <- m$conditions$estAngleMean
  zeroInflation <- m$conditions$zeroInflation
  oneInflation <- m$conditions$oneInflation
  DM <- m$conditions$DM
  DMind <- m$conditions$DMind
  
  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,m$conditions$bounds)
  
  mixtures <- m$conditions$mixtures
  
  fm <- NULL
  
  if(nbStates>1) {
    cat("Decoding state sequences for each imputation... \n")
    withCallingHandlers(im_states <- foreach(fm = im, i=seq_along(im), .combine = rbind) %dorng% {
      progBar(i,nsims)
      momentuHMM::viterbi(fm)
    },warning=muffleRNGwarning)
    if(nsims>1) states <- apply(im_states,2,function(x) which.max(hist(x,breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts))
    else states <- im_states
    cat("Decoding state probabilities for each imputation... \n")
    withCallingHandlers(im_stateProbs <- foreach(fm = im, i=seq_along(im)) %dorng% {
      progBar(i,nsims)
      momentuHMM::stateProbs(fm)
    },warning=muffleRNGwarning)
    if(mixtures>1){
      cat("Decoding mixture probabilities for each imputation... \n")
      withCallingHandlers(mixProbs <- foreach(fm = im, i=seq_along(im)) %dorng% {
        progBar(i,nsims)
        momentuHMM::mixtureProbs(fm)
      },warning=muffleRNGwarning)
    }
    #cat("DONE\n")
  } else states <- rep(1,nrow(data))
  if(ncores==1) doParallel::stopImplicitCluster()
  else future::plan(future::sequential)
  
  # pool estimates on working scale
  parms <- names(m$CIbeta)
  nparms <- length(parms)
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  parCount <- lapply(m$conditions$fullDM,ncol)#
  for(i in distnames[!unlist(lapply(m$conditions$circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
  }
  
  parmcols <- parCount
  parmcols$beta <- ncol(m$mle$beta)
  parmcols[["pi"]] <- mixtures-1
  parmcols$delta <- nbStates-1
  
  if(mixtures==1) piInd <- NULL
  else piInd <- 1
  
  parindex <- c(0,cumsum(c(unlist(parCount),length(m$mle$beta),piInd*ncol(m$covsPi)*(mixtures-1),ncol(m$covsDelta)*(nbStates-1)*mixtures)))
  names(parindex)[1:length(distnames)] <- distnames
  if(nbStates>1) {
    names(parindex)[length(distnames)+1] <- "beta"
    if(mixtures>1) names(parindex)[length(distnames)+2] <- "pi"
    names(parindex)[length(parindex)-1] <- "delta"
    if(!is.null(m$conditions$recharge)){
      parindex <- c(0,cumsum(c(unlist(parCount),length(m$mle$beta),piInd*ncol(m$covsPi)*(mixtures-1),ncol(m$covsDelta)*(nbStates-1)*mixtures,length(m$mle$g0),length(m$mle$theta))))
      if(mixtures>1) names(parindex)[1:(length(parindex)-1)] <- c(distnames,"beta","pi","delta","g0","theta")
      else names(parindex)[1:(length(parindex)-1)] <- c(distnames,"beta","delta","g0","theta")
      parmcols$g0 <- length(m$mle$g0)
      parmcols$theta <- length(m$mle$theta)
    }
  }
  parmcols <- unlist(parmcols[parms])
  
  miBeta <- mitools::MIcombine(results=betaCoeff,variances=betaVar)
  # account for betaCons
  if(nbStates>1){
    miBeta$variance[(parindex[["beta"]]+1:length(m$mle$beta))[duplicated(c(m$conditions$betaCons))],] <- 0
    miBeta$variance[,(parindex[["beta"]]+1:length(m$mle$beta))[duplicated(c(m$conditions$betaCons))]] <- 0
  }
  
  # multiple imputation results for working parameters
  twb <- lapply(im,function(x) x$mod$wpar)
  twb <- lapply(twb,function(x) {x[which(!is.finite(x))]<-NA; return(x)})
  if(length(m$conditions$optInd)) twvar <- lapply(im,function(x) x$mod$Sigma[-m$conditions$optInd,-m$conditions$optInd])
  else twvar <- lapply(im,function(x) x$mod$Sigma)
  twvar <- lapply(twvar,function(x) {x[which(!is.finite(x))]<-NA; diag(x)[which(!is.finite(sqrt(diag(x))))]<-NA; return(x)})
  miCombo <- mitools::MIcombine(results=twb,variances=twvar)
  
  for(parm in 1:nparms){
    
    parnames <- rownames(m$CIbeta[[parms[parm]]]$est)
    coeffs <- matrix(miBeta$coefficients[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    vars <- matrix(diag(miBeta$variance)[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    dfs <- matrix(miBeta$df[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    
    xbar[[parms[parm]]] <- matrix(NA,nrow=length(parnames),ncol=parmcols[parm])
    rownames(xbar[[parms[parm]]]) <- parnames
    MI_se[[parms[parm]]] <- lower[[parms[parm]]] <- upper[[parms[parm]]] <- xbar[[parms[parm]]]
    
    for(j in parnames){
      
      xbar[[parms[parm]]][j,] <- coeffs[j,]
      MI_se[[parms[parm]]][j,] <- sqrt(vars[j,])
      
      quantSup<-qt(1-(1-alpha)/2,df=dfs[j,])
      lower[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]-quantSup*MI_se[[parms[parm]]][j,]
      upper[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]+quantSup*MI_se[[parms[parm]]][j,]   
      
    }
  }
  
  Par <- list()
  Par$beta <- list()
  for(i in parms){
    Par$beta[[i]] <- mi_parm_list(xbar[[i]],MI_se[[i]],lower[[i]],upper[[i]],m$CIbeta[[i]]$est)
  }
  # fill in t.p.m. constraints based on betaCons
  if(nbStates>1) Par$beta$beta <- lapply(Par$beta$beta,function(x) matrix(x[c(m$conditions$betaCons)],dim(x),dimnames=list(rownames(x),colnames(x))))
  
  #average all numeric variables in imputed data
  mhdata<-m$data
  for(i in distnames){
    if(dist[[i]] %in% angledists) {
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,CircStats::circ.mean)
      class(mhdata[[i]]) <- c("angle",class(mhdata[[i]]))
    } else if(dist[[i]] %in% "pois"){
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,median)     
    } else if(dist[[i]] %in% mvndists){
      mhdata[[paste0(i,".x")]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[paste0(i,".x")]])),ncol=length(m$data[[paste0(i,".x")]]),byrow=TRUE),2,median)
      mhdata[[paste0(i,".y")]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[paste0(i,".y")]])),ncol=length(m$data[[paste0(i,".y")]]),byrow=TRUE),2,median)
      if(dist[[i]]=="mvnorm3" || dist[[i]]=="rw_mvnorm3")
        mhdata[[paste0(i,".z")]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[paste0(i,".z")]])),ncol=length(m$data[[paste0(i,".z")]]),byrow=TRUE),2,median)
    } else {
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,mean)
    }
  }
  for(j in names(m$data)[which(unlist(lapply(m$data,function(x) any(class(x) %in% meansListNoTime))) & !(names(m$data) %in% distnames))]){
    if(inherits(m$data[[j]],"angle")) {
      mhdata[[j]] <- apply(matrix(unlist(lapply(im,function(x) x$data[[j]])),ncol=length(m$data[[j]]),byrow=TRUE),2,CircStats::circ.mean)
      class(mhdata[[j]]) <- c("angle",class(mhdata[[j]]))
    } else mhdata[[j]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[j]])),ncol=length(m$data[[j]]),byrow=TRUE),2,mean)
  }
  mhrawCovs<-m$rawCovs
  if(length(mhrawCovs)){
    for(j in names(m$rawCovs)[which(unlist(lapply(m$rawCovs,function(x) any(class(x) %in% meansListNoTime))))]){
      mhrawCovs[[j]]<-apply(matrix(unlist(lapply(im,function(x) x$rawCovs[[j]])),ncol=length(m$rawCovs[[j]]),byrow=TRUE),2,mean)
    }
  }
  
  # identify covariates
  if(is.null(covs)){
    tempCovs <- mhdata[1,]
    for(j in names(mhdata)[which(unlist(lapply(mhdata,function(x) any(class(x) %in% meansList))))]){
      if(inherits(mhdata[[j]],"angle")) tempCovs[[j]] <- CircStats::circ.mean(mhdata[[j]][!is.na(mhdata[[j]])])
      else tempCovs[[j]]<-mean(mhdata[[j]],na.rm=TRUE)
    }
  } else {
    if(!is.data.frame(covs)) stop('covs must be a data frame')
    if(nrow(covs)>1) stop('covs must consist of a single row')
    if(!all(names(covs) %in% names(mhdata))) stop('invalid covs specified')
    if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(mhdata$ID))
    for(j in names(mhdata)[which(names(mhdata) %in% names(covs))]){
      if(inherits(mhdata[[j]],"factor")) covs[[j]] <- factor(covs[[j]],levels=levels(mhdata[[j]]))
      if(is.na(covs[[j]])) stop("check covs value for ",j)
    }    
    for(j in names(mhdata)[which(!(names(mhdata) %in% names(covs)))]){
      if(any(class(mhdata[[j]]) %in% meansList)) {
        if(inherits(mhdata[[j]],"angle")) covs[[j]] <- CircStats::circ.mean(mhdata[[j]][!is.na(mhdata[[j]])])
        else covs[[j]]<-mean(mhdata[[j]],na.rm=TRUE)
      } else covs[[j]] <- mhdata[[j]][1]
    }
    tempCovs <- covs[1,]
  }
  
  tmPar <- lapply(m$mle[distnames],function(x) c(t(x)))
  #parindex <- c(0,cumsum(unlist(parCount))[-length(m$conditions$fullDM)])
  #names(parindex) <- distnames
  for(i in distnames){
    if(!is.null(m$conditions$DM[[i]])){# & m$conditions$DMind[[i]]){
      tmPar[[i]] <- m$mod$estimate[parindex[[i]]+1:parCount[[i]]]
      if(!isFALSE(m$conditions$circularAngleMean[[i]])){
        names(tmPar[[i]]) <- unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]]))))
      } else names(tmPar[[i]])<-colnames(m$conditions$fullDM[[i]])
    } else if((dist[[i]] %in% angledists) & (!m$conditions$estAngleMean[[i]])){
      tmPar[[i]] <- tmPar[[i]][-(1:nbStates)]
    }
  }
  
  inputs <- checkInputs(nbStates,dist,tmPar,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$stateNames)
  p<-inputs$p
  splineInputs<-getSplineDM(distnames,inputs$DM,m,tempCovs)
  DMinputs<-getDM(splineInputs$covs,splineInputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,tmPar,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
  fullDM <- DMinputs$fullDM
  #DMinputs<-getDM(tempCovs,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,tmPar,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
  #fullDM<-DMinputs$fullDM
  
  # identify covariates
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,mhdata,covs=tempCovs,par=lapply(Par$beta,function(x) x$est))
  mhdata[colnames(reForm$newdata)] <- reForm$newdata
  attr(mhdata,'coords') <- attr(m$data,'coords')
  attr(mhdata,'coordLevel') <- attr(m$data,'coordLevel')
  recharge <- reForm$recharge
  hierRecharge <- reForm$hierRecharge
  newformula <- reForm$newformula
  tempCovs <- reForm$covs
  nbCovs <- reForm$nbCovs
  
  #miBeta <- mitools::MIcombine(results=lapply(im,function(x) x$mod$estimate),variances=lapply(im,function(x) x$mod$Sigma))
  
  ncmean <- get_ncmean(distnames,fullDM,m$conditions$circularAngleMean,nbStates)
  nc <- ncmean$nc
  meanind <- ncmean$meanind
  
  Par$real<-list()
  for(i in distnames){
    tmpParNames <- p$parNames[[i]]
    tmpParNames[which(p$parNames[[i]]=="kappa")] <- "concentration"
    
    DMind[[i]] <- FALSE
    par <- c(w2n(miBeta$coefficients,p$bounds,p$parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean[i],inputs$consensus[i],m$conditions$stationary,fullDM,DMind,1,inputs$dist[i],m$conditions$Bndind,nc,meanind,m$covsDelta,list(beta=matrix(rep(c(-Inf,Inf),length(m$mle$beta)),length(m$mle$beta),2,byrow=TRUE)),m$covsPi)[[i]])
    
    if(!(inputs$dist[[i]] %in% angledists) | (inputs$dist[[i]] %in% angledists & m$conditions$estAngleMean[[i]] & !m$conditions$Bndind[[i]])) {
      Par$real[[i]] <- get_CI(miBeta$coefficients,par,m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],NULL,inputs$dist[[i]])
    } else {
      if(!m$conditions$estAngleMean[[i]]){
        Par$real[[i]] <- get_CI(miBeta$coefficients,par[-(1:nbStates)],m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],NULL,inputs$dist[[i]])
        Par$real[[i]]$est <- matrix(c(rep(0,nbStates),Par$real[[i]]$est),ncol=nbStates,byrow=T)
        Par$real[[i]]$se <- matrix(c(rep(NA,nbStates),Par$real[[i]]$se),ncol=nbStates,byrow=T)
        Par$real[[i]]$lower <- matrix(c(rep(NA,nbStates),Par$real[[i]]$lower),ncol=nbStates,byrow=T)
        Par$real[[i]]$upper <- matrix(c(rep(NA,nbStates),Par$real[[i]]$upper),ncol=nbStates,byrow=T)  
        dimnames(Par$real[[i]]$est) <- dimnames(Par$real[[i]]$se) <- dimnames(Par$real[[i]]$lower) <- dimnames(Par$real[[i]]$upper) <- list(c("mean",tmpParNames),m$stateNames)
      } else {
        if(m$conditions$Bndind[[i]]){
          Par$real[[i]] <- CI_angle(miBeta$coefficients,par,m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],NULL,inputs$dist[[i]])
        }
      }
    }
  }
  
  quantSup<-qnorm(1-(1-alpha)/2)
  
  # pooled gamma estimates
  if(nbStates>1){
    gamInd <- (parindex[["beta"]]+1:((nbCovs+1)*nbStates*(nbStates-1)*mixtures))[unique(c(m$conditions$betaCons))]
    tmpSplineInputs<-getSplineFormula(newformula,mhdata,tempCovs)
    tempCovMat <- stats::model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
    
    est<-lower<-upper<-se<-matrix(NA,nbStates*mixtures,nbStates)
    
    for(mix in 1:mixtures){
      if(is.null(recharge)){
        wpar <- miBeta$coefficients[gamInd]
        est[(mix-1)*nbStates+1:nbStates,] <- get_gamma(wpar,tempCovMat,nbStates,1:nbStates,1:nbStates,m$conditions$betaRef,m$conditions$betaCons,mixture=mix)
        tmpSig <- miBeta$variance[gamInd,gamInd]
      } else {
        wpar <- c(miBeta$coefficients[gamInd],miBeta$coefficients[length(miBeta$coefficients)-reForm$nbRecovs:0])
        est[(mix-1)*nbStates+1:nbStates,] <- get_gamma_recharge(wpar,tmpSplineInputs$covs,tmpSplineInputs$formula,hierRecharge,nbStates,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,mixture=mix)
        tmpSig <- miBeta$variance[c(gamInd,length(miBeta$coefficients)-reForm$nbRecovs:0),c(gamInd,length(miBeta$coefficients)-reForm$nbRecovs:0)]
      }
      for(i in 1:nbStates){
        for(j in 1:nbStates){
          if(is.null(recharge)){
            dN<-numDeriv::grad(get_gamma,wpar,covs=tempCovMat,nbStates=nbStates,i=i,j=j,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,mixture=mix)
          } else {
            dN<-numDeriv::grad(get_gamma_recharge,wpar,covs=tmpSplineInputs$covs,formula=tmpSplineInputs$formula,hierRecharge=hierRecharge,nbStates=nbStates,i=i,j=j,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,mixture=mix)
          }  
          se[(mix-1)*nbStates+i,j]<-suppressWarnings(sqrt(dN%*%tmpSig%*%dN))
          lower[(mix-1)*nbStates+i,j]<-1/(1+exp(-(log(est[(mix-1)*nbStates+i,j]/(1-est[(mix-1)*nbStates+i,j]))-quantSup*(1/(est[(mix-1)*nbStates+i,j]-est[(mix-1)*nbStates+i,j]^2))*se[(mix-1)*nbStates+i,j])))#est[(mix-1)*nbStates+i,j]-quantSup*se[(mix-1)*nbStates+i,j]
          upper[(mix-1)*nbStates+i,j]<-1/(1+exp(-(log(est[(mix-1)*nbStates+i,j]/(1-est[(mix-1)*nbStates+i,j]))+quantSup*(1/(est[(mix-1)*nbStates+i,j]-est[(mix-1)*nbStates+i,j]^2))*se[(mix-1)*nbStates+i,j])))#est[(mix-1)*nbStates+i,j]+quantSup*se[(mix-1)*nbStates+i,j]
        }
      }
    }
    Par$real$gamma <- list(est=est,se=se,lower=lower,upper=upper)
    dimnames(Par$real$gamma$est) <- dimnames(Par$real$gamma$se) <- dimnames(Par$real$gamma$lower) <- dimnames(Par$real$gamma$upper) <- list(rep(m$stateNames,mixtures),m$stateNames)
    if(mixtures>1) dimnames(Par$real$gamma$est) <- dimnames(Par$real$gamma$se) <- dimnames(Par$real$gamma$lower) <- dimnames(Par$real$gamma$upper) <- list(paste0(rep(m$stateNames,mixtures),"_mix",rep(1:mixtures,each=nbStates)),m$stateNames)
  }
  
  # pooled pi estimates
  if(mixtures>1 & nbStates>1){
    piInd <- parindex[["beta"]]+((nbCovs+1)*nbStates*(nbStates-1)*mixtures)+1:(ncol(m$covsPi)*(mixtures-1))
    pie <- matrix(miBeta$coefficients[piInd],nrow=ncol(m$covsPi),ncol=mixtures-1)
    est<-lower<-upper<-se<-matrix(NA,nrow=nrow(m$covsPi),ncol=mixtures)
    for(j in 1:nrow(m$covsPi)){
      est[j,] <- get_delta(pie,m$covsPi[j,,drop=FALSE],i=1:mixtures)
      for(i in 1:mixtures){
        dN<-numDeriv::grad(get_delta,pie,covsDelta=m$covsPi[j,,drop=FALSE],i=i)
        se[j,i]<-suppressWarnings(sqrt(dN%*%miBeta$variance[piInd,piInd]%*%dN))
        lower[j,i] <- probCI(est[j,i],se[j,i],quantSup,bound="lower")
        upper[j,i] <- probCI(est[j,i],se[j,i],quantSup,bound="upper")
      }
    }
    Par$real[["pi"]] <- list(est=est,se=se,lower=lower,upper=upper)
    colnames(Par$real[["pi"]]$est) <- colnames(Par$real[["pi"]]$se) <- colnames(Par$real[["pi"]]$lower) <- colnames(Par$real[["pi"]]$upper) <- paste0("mix",1:mixtures)
    rownames(Par$real[["pi"]]$est) <- rownames(Par$real[["pi"]]$se) <- rownames(Par$real[["pi"]]$lower) <- rownames(Par$real[["pi"]]$upper) <- paste0("ID:",unique(m$data$ID))
  }
  
  # pooled delta estimates
  if(!m$conditions$stationary & nbStates>1){
    nbCovsDelta <- ncol(m$covsDelta)-1
    foo <- length(miBeta$coefficients)-ifelse(reForm$nbRecovs,(reForm$nbRecovs+1)+(reForm$nbG0covs+1),0)-(nbCovsDelta+1)*(nbStates-1)*mixtures
    deltInd <- foo+1:((nbCovsDelta+1)*(nbStates-1)*mixtures)
    delta <- matrix(miBeta$coefficients[deltInd],nrow=(nbCovsDelta+1)*mixtures,ncol=nbStates-1)
    est<-lower<-upper<-se<-matrix(NA,nrow=nrow(m$covsDelta)*mixtures,ncol=nbStates)
    for(mix in 1:mixtures){
      for(j in 1:nrow(m$covsDelta)){
        est[(mix-1)*nrow(m$covsDelta)+j,] <- get_delta(delta,m$covsDelta[j,,drop=FALSE],1:nbStates,mixture=mix)
        for(i in 1:nbStates){
          dN<-numDeriv::grad(get_delta,delta,covsDelta=m$covsDelta[j,,drop=FALSE],i=i,mixture=mix)
          se[(mix-1)*nrow(m$covsDelta)+j,i]<-suppressWarnings(sqrt(dN%*%miBeta$variance[deltInd,deltInd]%*%dN))
          lower[(mix-1)*nrow(m$covsDelta)+j,i] <- probCI(est[(mix-1)*nrow(m$covsDelta)+j,i],se[(mix-1)*nrow(m$covsDelta)+j,i],quantSup,bound="lower")
          upper[(mix-1)*nrow(m$covsDelta)+j,i] <- probCI(est[(mix-1)*nrow(m$covsDelta)+j,i],se[(mix-1)*nrow(m$covsDelta)+j,i],quantSup,bound="upper")
        }
      }
    }
  } else {
    if(nbStates>1){
      covs<-stats::model.matrix(newformula,tempCovs)
      statFun<-function(beta,nbStates,covs,i,mixture=1){
        gamma <- trMatrix_rcpp(nbStates,beta[(mixture-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],covs,m$conditions$betaRef)[,,1]
        tryCatch(solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))[i],error = function(e) {
          "A problem occurred in the calculation of the stationary distribution."})
      }
      est <- lower <- upper <- se <- matrix(NA,nbAnimals*mixtures,nbStates)
      for(mix in 1:mixtures){
        delta <- statFun(matrix(miBeta$coefficients[gamInd],nrow=(nbCovs+1)*mixtures),nbStates,covs,1:nbStates,mixture=mix)
        est[nbAnimals*(mix-1)+1:nbAnimals,] <- matrix(delta,nrow=nbAnimals,ncol=nbStates,byrow=TRUE)
        for(k in 1:nbStates){
          dN<-numDeriv::grad(statFun,matrix(miBeta$coefficients[gamInd],nrow=(nbCovs+1)*mixtures),nbStates=nbStates,covs=covs,i=k,mixture=mix)
          se[nbAnimals*(mix-1)+1:nbAnimals,k]<-suppressWarnings(sqrt(dN%*%miBeta$variance[gamInd,gamInd]%*%dN))
          lower[nbAnimals*(mix-1)+1:nbAnimals,k] <- probCI(est[nbAnimals*(mix-1)+1:nbAnimals,k],se[nbAnimals*(mix-1)+1:nbAnimals,k],quantSup,bound="lower")
          upper[nbAnimals*(mix-1)+1:nbAnimals,k] <- probCI(est[nbAnimals*(mix-1)+1:nbAnimals,k],se[nbAnimals*(mix-1)+1:nbAnimals,k],quantSup,bound="upper")
        }
      }
    } else {
      est <- matrix(1,nrow(m$covsDelta)*mixtures)
      lower <- upper <- se <- matrix(NA,nrow(m$covsDelta))
    }
  }
  Par$real$delta <- list(est=est,se=se,lower=lower,upper=upper)
  colnames(Par$real$delta$est) <- colnames(Par$real$delta$se) <- colnames(Par$real$delta$lower) <- colnames(Par$real$delta$upper) <- m$stateNames
  rownames(Par$real$delta$est) <- rownames(Par$real$delta$se) <- rownames(Par$real$delta$lower) <- rownames(Par$real$delta$upper) <- paste0("ID:",rep(unique(m$data$ID),mixtures))
  if(mixtures>1) rownames(Par$real$delta$est) <- rownames(Par$real$delta$se) <- rownames(Par$real$delta$lower) <- rownames(Par$real$delta$upper) <- paste0("ID:",rep(unique(m$data$ID),mixtures),"_mix",rep(1:mixtures,each=nbAnimals))
  
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  
  if(nbStates>1){
    xmat[["stateProbs"]] <- array(unlist(im_stateProbs),c(nrow(data),nbStates,nsims))
    xvar[["stateProbs"]] <- array(0,c(nrow(data),nbStates,nsims)) # don't have se's; might be a way to get these but probably quite complicated
    n <- apply(!(is.na(xmat[["stateProbs"]])+is.na(xvar[["stateProbs"]])),1:2,sum)
    
    if(any(n<2)) warning("need at least 2 simulations with valid point and variance estimates for stateProbs")
    
    xbar[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    B_m[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , var,na.rm=TRUE)
    
    W_m[["stateProbs"]] <- apply( xvar[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    MI_se[["stateProbs"]] <- sqrt(W_m[["stateProbs"]] + (n+1)/n * B_m[["stateProbs"]])
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["stateProbs"]]/B_m[["stateProbs"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    lower[["stateProbs"]] <- suppressWarnings(probCI(xbar[["stateProbs"]],MI_se[["stateProbs"]],quantSup,"lower"))
    upper[["stateProbs"]] <- suppressWarnings(probCI(xbar[["stateProbs"]],MI_se[["stateProbs"]],quantSup,"upper"))
    
    xmat[["timeInStates"]] <- t(apply(im_states,1,function(x) {counts<-hist(x,breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts;counts/sum(counts)}))
    xvar[["timeInStates"]] <- matrix(0 , ncol=nbStates, nrow=nsims, byrow=TRUE) # don't have se's; might be a way to get these but probably quite complicated
    n <- apply(!(is.na(xmat[["timeInStates"]])+is.na(xvar[["timeInStates"]])),2,sum)
    
    if(any(n<2)) warning("need at least 2 simulations with valid point and variance estimates for timeInStates")
    
    xbar[["timeInStates"]] <- apply(xmat[["timeInStates"]],2,mean,na.rm=TRUE)
    B_m[["timeInStates"]] <- apply(xmat[["timeInStates"]],2,var,na.rm=TRUE)
    
    W_m[["timeInStates"]] <- apply(xvar[["timeInStates"]],2,mean,na.rm=TRUE)
    MI_se[["timeInStates"]] <- sqrt(W_m[["timeInStates"]] + (n+1)/n * B_m[["timeInStates"]])
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["timeInStates"]]/B_m[["timeInStates"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    lower[["timeInStates"]] <- probCI(xbar[["timeInStates"]],MI_se[["timeInStates"]],quantSup,"lower")
    upper[["timeInStates"]] <- probCI(xbar[["timeInStates"]],MI_se[["timeInStates"]],quantSup,"upper")
    
    if(mixtures>1){
      xmat[["mixtureProbs"]] <- array(unlist(mixProbs),c(nbAnimals,mixtures,nsims))
      xvar[["mixtureProbs"]] <- array(0,c(nbAnimals,mixtures,nsims)) # don't have se's; might be a way to get these but probably quite complicated
      n <- apply(!(is.na(xmat[["mixtureProbs"]])+is.na(xvar[["mixtureProbs"]])),1:2,sum)
      
      if(any(n<2)) warning("need at least 2 simulations with valid point and variance estimates for mixtureProbs")
      
      xbar[["mixtureProbs"]] <-   apply( xmat[["mixtureProbs"]] , 1:2 , mean,na.rm=TRUE)
      B_m[["mixtureProbs"]] <-   apply( xmat[["mixtureProbs"]] , 1:2 , var,na.rm=TRUE)
      
      W_m[["mixtureProbs"]] <- apply( xvar[["mixtureProbs"]] , 1:2 , mean,na.rm=TRUE)
      MI_se[["mixtureProbs"]] <- sqrt(W_m[["mixtureProbs"]] + (n+1)/n * B_m[["mixtureProbs"]])
      
      dfs<-(n-1)*(1+1/(n+1)*W_m[["mixtureProbs"]]/B_m[["mixtureProbs"]])^2
      quantSup<-qt(1-(1-alpha)/2,df=dfs)
      
      lower[["mixtureProbs"]] <- suppressWarnings(probCI(xbar[["mixtureProbs"]],MI_se[["mixtureProbs"]],quantSup,"lower"))
      upper[["mixtureProbs"]] <- suppressWarnings(probCI(xbar[["mixtureProbs"]],MI_se[["mixtureProbs"]],quantSup,"upper"))
    }
  }
  
  if(nbStates>1) {
    Par$timeInStates <- list(est=xbar$timeInStates,se=MI_se$timeInStates,lower=lower$timeInStates,upper=upper$timeInStates)
    Par$timeInStates <- lapply(Par$timeInStates,function(x){ names(x) = m$stateNames;x})
    
    Par$states <- states
    
    Par$stateProbs <- list(est=xbar$stateProbs,se=MI_se$stateProbs,lower=lower$stateProbs,upper=upper$stateProbs)
    Par$stateProbs <- lapply(Par$stateProbs,function(x) {rownames(x) = data$ID;x})
    Par$stateProbs <- lapply(Par$stateProbs,function(x) {colnames(x) = m$stateNames;x})
    
    if(mixtures>1){
      Par$mixtureProbs <- list(est=xbar$mixtureProbs,se=MI_se$mixtureProbs,lower=lower$mixtureProbs,upper=upper$mixtureProbs)
      Par$mixtureProbs <- lapply(Par$mixtureProbs,function(x) {rownames(x)=paste0("ID:",unique(data$ID));x})
      Par$mixtureProbs <- lapply(Par$mixtureProbs,function(x) {colnames(x)=paste0("mix",1:mixtures);x})
    }
  }
  
  if(inherits(im[[1]],"hierarchical")){
    tmp<-lapply(Par$stateProbs,function(x) hierStateProbs(im[[1]],x))
    Par$hierStateProbs <- list()
    for(j in names(tmp$est)){
      Par$hierStateProbs[[j]] <- list()
      for(jj in names(tmp)){
        Par$hierStateProbs[[j]][[jj]] <- tmp[[jj]][[j]]
      }
    }
  }
  
  mh <- im[[1]]
  attr(mh,"class") <- NULL
  mh$mle <- NULL
  mh$mod <- NULL
  mh$CIreal <- NULL
  mh$CIbeta <- NULL
  if(any(ident)) mh$conditions$fullDM <- fullDM
  
  mh$data<-mhdata[!(colnames(mhdata) %in% colnames(reForm$newdata))]
  mh$rawCovs<-mhrawCovs
  
  # get fixPar$delta in working scale format so expandPar works correctly in post-analysis
  if(!is.momentuHierHMM(im[[1]]) & !mh$conditions$stationary & nbStates>1) {
    if(any(!is.na(mh$conditions$fixPar$delta))){
      tmp <- which(!is.na(mh$conditions$fixPar$delta))
      if(!nbCovsDelta){
        delta0 <- mh$conditions$fixPar$delta
        delta0 <- matrix(delta0,mixtures,nbStates)
        delta0 <- matrix(apply(delta0,1,function(x) log(x[-1]/x[1])),(nbCovsDelta+1)*mixtures,nbStates-1)
        mh$conditions$fixPar$delta <- as.vector(delta0)
      }
    } else mh$conditions$fixPar$delta <- rep(NA,length(deltInd))
  }
  
  coordNames <- attr(m$data,"coords")
  
  mvnorm2Ind <- 1
  if(!is.null(m$conditions$mvnCoords)){
    coordNames <- c("x","y")
    if(m$conditions$dist[[m$conditions$mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) mvnorm2Ind <- 0#coordNames <- c("x","y","z")
    coordNames <- paste0(m$conditions$mvnCoords,".",coordNames)
  } else if(is.null(coordNames)) coordNames <- c("x","y")
  
  errorEllipse<-NULL
  if(all(coordNames %in% names(mh$data)) & mvnorm2Ind){
    checkerrs <- lapply(im,function(x) x$data[match(coordNames,names(x$data))])
    ident <- !unlist(lapply(checkerrs,function(x) isTRUE(all.equal(x,checkerrs[[1]]))))
    if(any(ident)){
      # calculate location alpha% error ellipses
      if (!requireNamespace("car", quietly = TRUE)) {
        warning("Package \"car\" needed for calculating error ellipses. Please install it.",
             call. = FALSE)
      } else {
        if(ncores>1){
          future::plan(future::multisession, workers = ncores)
        } else { 
          doParallel::registerDoParallel(cores=ncores)
        }
        cat("Calculating location",paste0(alpha*100,"%"),"error ellipses... ")
        tmpx<-matrix(unlist(lapply(im,function(x) x$data[[coordNames[1]]])),nrow(mh$data))
        tmpy<-matrix(unlist(lapply(im,function(x) x$data[[coordNames[2]]])),nrow(mh$data))
        withCallingHandlers(errorEllipse<-foreach(i = 1:nrow(mh$data)) %dorng% {
          tmp <- cbind(tmpx[i,],tmpy[i,])
          if(length(unique(tmp[,1]))>1 | length(unique(tmp[,2]))>1)
            ellip <- car::dataEllipse(tmp,levels=alpha,draw=FALSE,segments=100)
          else ellip <- matrix(tmp[1,],101,2,byrow=TRUE)
        },warning=muffleRNGwarning)
        if(ncores==1) doParallel::stopImplicitCluster()
        else future::plan(future::sequential)
        cat("DONE\n")
      }
    }
  }
  
  mh$errorEllipse <- errorEllipse
  mh$Par <- Par
  mh$MIcombine <- miCombo
  
  mh <- miSum(mh)
  if(is.momentuHierHMM(im[[1]])) class(mh) <- append(class(mh),"hierarchical")
  
  if(inherits(mh,"hierarchical")){
    inputHierHMM <- formatHierHMM(mh$data,mh$conditions$hierStates,mh$conditions$hierDist,hierBeta=NULL,hierDelta=NULL,mh$conditions$hierFormula,mh$conditions$hierFormulaDelta,mh$conditions$mixtures)
    hier <- mapHier(list(beta=mh$Par$beta$beta$est,g0=mh$Par$beta$g0$est,theta=mh$Par$beta$theta$est),mh$Par$beta[["pi"]]$est,mh$Par$beta$delta$est,mh$conditions$hierBeta,mh$conditions$hierDelta,inputHierHMM$hFixPar,inputHierHMM$hBetaCons,inputHierHMM$hDeltaCons,mh$conditions$hierStates,inputHierHMM$newformula,mh$conditions$formulaDelta,inputHierHMM$data,mh$conditions$mixtures,inputHierHMM$recharge)
    mh$conditions$hierBeta <- hier$hierBeta
    mh$conditions$hierDelta <- hier$hierDelta
    
    mh$Par$real <- CIreal.hierarchical(mh)
  }
  
  attr(mh$data,"coords") <- coordNames
  
  return(mh)
}

mi_parm_list<-function(est,se,lower,upper,m){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  dimnames(Par$est) <- dimnames(Par$se) <- dimnames(Par$lower) <- dimnames(Par$upper) <- list(rownames(m),colnames(m))
  Par
}

#' @importFrom stats plogis qlogis
probCI<-function(x,se,z,bound="lower"){
  if(bound=="lower")
    ci<-stats::plogis(stats::qlogis(x)-z*(1/(x-x^2))*se) 
  else if(bound=="upper")
    ci<-stats::plogis(stats::qlogis(x)+z*(1/(x-x^2))*se) 
  ci
}

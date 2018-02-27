#'
#' Calculate pooled parameter estimates and states across multiple imputations
#' 
#' @param HMMfits List comprised of \code{\link{momentuHMM}} objects
#' @param alpha Significance level for calculating confidence intervals of pooled estimates (including location error ellipses). Default: 0.95.
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing).
#' @param covs Data frame consisting of a single row indicating the covariate values to be used in the calculation of pooled natural parameters. 
#' For any covariates that are not specified using \code{covs}, the means of the covariate(s) across the imputations are used 
#' (unless the covariate is a factor, in which case the first factor in the data is used). By default, no covariates are specified.
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
#' 
#' @details
#' Pooled estimates, standard errors, and confidence intervals are calculated using standard multiple imputation formulas. Working scale parameters are pooled
#' using \code{\link[mitools]{MIcombine}} and t-distributed confidence intervals. Natural scale parameters and normally-distributed confidence intervals are calculated by transforming the pooled working scale parameters 
#' and, if applicable, are based on covariate means across all imputations (and/or values specified in \code{covs}).
#' 
#' Note that pooled estimates for \code{timeInStates} and \code{stateProbs} do not include within-model uncertainty and are based entirely on across-model variability.
#' 
#' @examples
#' \dontrun{
#' # Extract data and crawl inputs from miExample
#' obsData <- miExample$obsData
#' inits <- miExample$inits
#' err.model <- miExample$err.model
#' 
#' # Fit crawl to obsData
#' crwOut <- crawlWrap(obsData,theta=c(4,0),fixPar=c(1,1,NA,NA),
#'                     initial.state=inits,err.model=err.model)
#'                     
#' # Fit four imputations
#' bPar <- miExample$bPar
#' HMMfits <- MIfitHMM(crwOut,nSims=4,poolEstimates=FALSE,
#'                    nbStates=2,dist=list(step="gamma",angle="vm"),
#'                    Par0=bPar$Par,beta0=bPar$beta,delta0=bPar$delta,
#'                    formula=~cov1+cos(cov2),
#'                    estAngleMean=list(angle=TRUE),
#'                    covNames=c("cov1","cov2"))
#'                    
#' # Pool estimates
#' miSum <- MIpool(HMMfits)
#' print(miSum)
#' }
#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats median var qt
#' @importFrom boot logit inv.logit
#' @importFrom CircStats circ.mean
#' @importFrom car dataEllipse
#' @importFrom mitools MIcombine
#' @importFrom MASS ginv
MIpool<-function(HMMfits,alpha=0.95,ncores=1,covs=NULL){
  
  im <- HMMfits
  goodIndex <- 1:length(im)
  simind <- which((unlist(lapply(im,is.momentuHMM))))
  nsims <- length(simind)
  if(nsims<1) stop("'HMMfits' must be a list comprised of momentuHMM objects")
  
  checkmove <- which(!(unlist(lapply(im,is.momentuHMM))))
  if(length(checkmove)) {
    im[checkmove]<-NULL
    warning("The following imputations are not momentuHMM objects and will be ignored: ",paste(checkmove,collapse=", "))
    goodIndex <- goodIndex[-checkmove]
  }
  checksims <- lapply(im,function(x) x[match("conditions",names(x))])
  ident <- !unlist(lapply(checksims,function(x) isTRUE(all.equal(x,checksims[[1]]))))
  if(any(ident)){
    # check that only differences are in the design matrix covariate values
    checksims2 <- lapply(checksims, function(x) x$conditions[-match("fullDM",names(x$conditions))])
    ident2 <- !unlist(lapply(checksims2,function(x) isTRUE(all.equal(x,checksims2[[1]]))))
    if(any(ident2)) stop("Model conditions for each imputation must be identical. Imputations that do not match the first: ",paste(which(ident),collapse=", "))
  }
  
  if(any(unlist(lapply(im,function(x) is.null(x$mod$hessian))))) stop("Estimates cannot be pooled unless Hessian is calculated. Hessian is missing for imputations ",paste0(which(unlist(lapply(im,function(x) is.null(x$mod$hessian)))),collapse=", "))
  
  tmpDet <- which(unlist(lapply(im,function(x) det(x$mod$hessian)))==0)
  if(length(tmpDet)){
    warning("Hessian is singular for HMM fit(s): ",paste0(goodIndex[tmpDet],collapse=", "))
  }
  
  tmpVar <- which(unlist(lapply(im,function(x) any(class(tryCatch(ginv(x$mod$hessian),error=function(e) e)) %in% "error"))))
  if(length(tmpVar)){
    warning("ginv of the hessian failed for HMM fit(s): ",paste0(goodIndex[tmpVar],collapse=", "))
    im[tmpVar] <- NULL
    nsims <- length(im)
    goodIndex <- goodIndex[-tmpVar]
  }
  
  im <- lapply(im,delta_bc)
  m <- im[[1]]
  
  tempcons<-rep(1,length(m$mod$estimate))
  tempworkcons<-rep(0,length(m$mod$estimate))
  tempcons[1:length(unlist(m$conditions$cons))]<-unlist(m$conditions$cons)
  tempworkcons[1:length(unlist(m$conditions$workcons))]<-unlist(m$conditions$workcons)
  
  wBounds <- cbind(unlist(lapply(m$conditions$workBounds,function(x) x[,1])),unlist(lapply(m$conditions$workBounds,function(x) x[,2])))
  
  # check for finite coefficients and standard errors
  betaVar <- lapply(im,function(x) get_gradwb(x$mod$estimate,wBounds,tempcons)%*%ginv(x$mod$hessian)%*%t(get_gradwb(x$mod$estimate,wBounds,tempcons)))
  betaCoeff <- lapply(im,function(x) w2wn(x$mod$estimate^tempcons+tempworkcons,wBounds))
  tmpVar1 <- which(unlist(lapply(betaCoeff,function(x) any(!is.finite(x)))))
  if(length(tmpVar1)){
    warning("working parameter estimates are not finite for HMM fits ",paste0(goodIndex[tmpVar1],collapse=", ")," and will not be included in pooling")
    im[tmpVar1] <- NULL
    m <- im[[1]]
    nsims <- length(im)
    betaVar <- lapply(im,function(x) get_gradwb(x$mod$estimate,wBounds,tempcons)%*%ginv(x$mod$hessian)%*%t(get_gradwb(x$mod$estimate,wBounds,tempcons)))
    betaCoeff <- lapply(im,function(x) w2wn(x$mod$estimate^tempcons+tempworkcons,wBounds))
    goodIndex <- goodIndex[-tmpVar1]
  }
  tmpVar2 <- which(unlist(lapply(betaVar,function(x) any(!is.finite(x)))))
  if(length(tmpVar2)){
    warning("working parameter standard errors are not finite for HMM fits ",paste0(goodIndex[tmpVar2],collapse=", ")," and will not be included in pooling")
    im[tmpVar2] <- NULL
    m <- im[[1]]
    nsims <- length(im)
    betaCoeff <- lapply(im,function(x) w2wn(x$mod$estimate^tempcons+tempworkcons,wBounds))
    betaVar <- lapply(im,function(x) get_gradwb(x$mod$estimate,wBounds,tempcons)%*%ginv(x$mod$hessian)%*%t(get_gradwb(x$mod$estimate,wBounds,tempcons)))
    goodIndex <- goodIndex[-tmpVar2]
  }
  
  data <- m$data
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  estAngleMean <- m$conditions$estAngleMean
  zeroInflation <- m$conditions$zeroInflation
  oneInflation <- m$conditions$oneInflation
  DM <- m$conditions$DM
  DMind <- m$conditions$DMind
  
  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,m$conditions$bounds)
  
  if(nbStates>1) {
    cat("Decoding state sequences and probabilities for each imputation... ")
    registerDoParallel(cores=ncores)
    im_states <- foreach(i = 1:nsims, .combine = rbind) %dopar% {momentuHMM::viterbi(im[[i]])}
    stopImplicitCluster()
    states <- apply(im_states,2,function(x) which.max(hist(x,breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts))
    registerDoParallel(cores=ncores)
    im_stateProbs <- foreach(i = 1:nsims) %dopar% {momentuHMM::stateProbs(im[[i]])}
    stopImplicitCluster()
    cat("DONE\n")
  } else states <- rep(1,nrow(data))
  
  # pool estimates on working scale
  parms <- names(m$CIbeta)
  nparms <- length(parms)
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  parCount <- lapply(m$conditions$fullDM,ncol)#
  for(i in distnames[unlist(m$conditions$circularAngleMean)]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
  }
  parmcols <- parCount
  parmcols$beta <- ncol(m$mle$beta)
  parmcols$delta <- nbStates-1
  parmcols <- unlist(parmcols[parms])
  
  parindex <- c(0,cumsum(c(unlist(parCount),length(m$mle$beta),ncol(m$covsDelta)*(nbStates-1))))
  names(parindex)[1:length(distnames)] <- distnames
  if(nbStates>1) {
    names(parindex)[length(distnames)+1] <- "beta"
    names(parindex)[length(parindex)-1] <- "delta"
  }
  
  miBeta <- mitools::MIcombine(results=betaCoeff,variances=betaVar)
  
  for(parm in 1:nparms){
    
    parnames <- rownames(m$CIbeta[[parms[parm]]]$est)
    #if(parms[parm] %in% distnames){
    #  coeffs <- matrix(miBeta$coefficients[(parindex[parm]+1):parindex[parm+1]]^m$conditions$cons[[parms[parm]]]+m$conditions$workcons[[parms[parm]]],nrow=length(parnames),dimnames=list(parnames))
    #  vars <- matrix(((m$conditions$cons[[parms[parm]]]*(miBeta$coefficients[(parindex[parm]+1):parindex[parm+1]]^(m$conditions$cons[[parms[parm]]]-1)))^2)*(diag(miBeta$variance)[(parindex[parm]+1):parindex[parm+1]]),nrow=length(parnames),dimnames=list(parnames))
    #} else {
    coeffs <- matrix(miBeta$coefficients[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    vars <- matrix(diag(miBeta$variance)[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    #}
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
  
  #average all numeric variables in imputed data
  mhdata<-m$data
  for(i in distnames){
    if(dist[[i]] %in% angledists) {
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,CircStats::circ.mean)
    } else if(dist[[i]] %in% "pois"){
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,median)     
    } else {
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,mean)
    }
  }
  for(j in names(m$data)[which(unlist(lapply(m$data,function(x) any(class(x) %in% meansListNoTime))) & !(names(m$data) %in% distnames))]){
    mhdata[[j]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[j]])),ncol=length(m$data[[j]]),byrow=TRUE),2,mean)
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
  parindex <- c(0,cumsum(unlist(parCount))[-length(m$conditions$fullDM)])
  names(parindex) <- distnames
  for(i in distnames){
    if(!is.null(m$conditions$DM[[i]])){# & m$conditions$DMind[[i]]){
      tmPar[[i]] <- m$mod$estimate[parindex[[i]]+1:parCount[[i]]]
      if(m$conditions$circularAngleMean[[i]]){
        names(tmPar[[i]]) <- unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]]))))
      } else names(tmPar[[i]])<-colnames(m$conditions$fullDM[[i]])
    } else if((dist[[i]] %in% angledists) & (!m$conditions$estAngleMean[[i]])){
      tmPar[[i]] <- tmPar[[i]][-(1:nbStates)]
    }
  }
  
  inputs <- checkInputs(nbStates,dist,tmPar,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,NULL,NULL,m$stateNames)
  p<-inputs$p
  splineInputs<-getSplineDM(distnames,inputs$DM,m,tempCovs)
  DMinputs<-getDM(splineInputs$covs,splineInputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,tmPar,NULL,NULL,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
  fullDM <- DMinputs$fullDM
  #DMinputs<-getDM(tempCovs,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,tmPar,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
  #fullDM<-DMinputs$fullDM
  
  formula<-m$conditions$formula
  stateForms<- terms(formula, specials = paste0("state",1:nbStates))
  newformula<-formula
  if(nbStates>1){
    if(length(unlist(attr(stateForms,"specials")))){
      newForm<-attr(stateForms,"term.labels")[-unlist(attr(stateForms,"specials"))]
      for(i in 1:nbStates){
        if(!is.null(attr(stateForms,"specials")[[paste0("state",i)]])){
          for(j in 1:(nbStates-1)){
            newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("state",i)]]]))
          }
        }
      }
      newformula<-as.formula(paste("~",paste(newForm,collapse="+")))
    }
    formulaStates<-stateFormulas(newformula,nbStates*(nbStates-1),spec="betaCol")
    if(length(unlist(attr(terms(newformula, specials = c(paste0("betaCol",1:(nbStates*(nbStates-1))),"cosinor")),"specials")))){
      allTerms<-unlist(lapply(formulaStates,function(x) attr(terms(x),"term.labels")))
      newformula<-as.formula(paste("~",paste(allTerms,collapse="+")))
      formterms<-attr(terms.formula(newformula),"term.labels")
    } else {
      formterms<-attr(terms.formula(newformula),"term.labels")
      newformula<-formula
    }
  }
  nbCovs <- ncol(model.matrix(newformula,m$data))-1 # substract intercept column
  
  #miBeta <- mitools::MIcombine(results=lapply(im,function(x) x$mod$estimate),variances=lapply(im,function(x) ginv(x$mod$hessian)))
  
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) {
      meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  
  Par$real<-list()
  for(i in distnames){
    tmpParNames <- p$parNames[[i]]
    tmpParNames[which(p$parNames[[i]]=="kappa")] <- "concentration"
    
    DMind[[i]] <- FALSE
    par <- c(w2n(miBeta$coefficients,p$bounds,p$parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean[i],inputs$consensus[i],m$conditions$stationary,DMinputs$cons,fullDM,DMind,DMinputs$workcons,1,inputs$dist[i],m$conditions$Bndind,nc,meanind,m$covsDelta,NULL)[[i]])

    if(!(inputs$dist[[i]] %in% angledists) | (inputs$dist[[i]] %in% angledists & m$conditions$estAngleMean[[i]] & !m$conditions$Bndind[[i]])) {
      Par$real[[i]] <- get_CI(miBeta$coefficients,par,m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],DMinputs$cons[[i]],DMinputs$workcons[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],NULL)
    } else {
      if(!m$conditions$estAngleMean[[i]]){
        Par$real[[i]] <- get_CI(miBeta$coefficients,par[-(1:nbStates)],m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],DMinputs$cons[[i]],DMinputs$workcons[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],NULL)
        Par$real[[i]]$est <- matrix(c(rep(0,nbStates),Par$real[[i]]$est),ncol=nbStates,byrow=T)
        Par$real[[i]]$se <- matrix(c(rep(NA,nbStates),Par$real[[i]]$se),ncol=nbStates,byrow=T)
        Par$real[[i]]$lower <- matrix(c(rep(NA,nbStates),Par$real[[i]]$lower),ncol=nbStates,byrow=T)
        Par$real[[i]]$upper <- matrix(c(rep(NA,nbStates),Par$real[[i]]$upper),ncol=nbStates,byrow=T)  
        dimnames(Par$real[[i]]$est) <- dimnames(Par$real[[i]]$se) <- dimnames(Par$real[[i]]$lower) <- dimnames(Par$real[[i]]$upper) <- list(c("mean",tmpParNames),m$stateNames)
      } else {
        if(m$conditions$Bndind[[i]]){
          Par$real[[i]] <- CI_angle(miBeta$coefficients,par,m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],DMinputs$cons[[i]],DMinputs$workcons[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],NULL)
        }
      }
    }
  }
  
  quantSup<-qnorm(1-(1-alpha)/2)
    
  # pooled gamma estimates
  if(nbStates>1){
    gamInd<-(length(miBeta$coefficients)-(nbCovs+1)*nbStates*(nbStates-1)+1):(length(miBeta$coefficients))-ncol(m$covsDelta)*(nbStates-1)*(1-m$conditions$stationary)
    tmpSplineInputs<-getSplineFormula(newformula,mhdata,tempCovs)
    tempCovMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
    est <- get_gamma(matrix(miBeta$coefficients[gamInd],nrow=nbCovs+1),tempCovMat,nbStates,1:nbStates,1:nbStates)
    lower<-upper<-se<-matrix(0,nrow(est),ncol(est))
    for(i in 1:nrow(est)){
      for(j in 1:ncol(est)){
        dN<-numDeriv::grad(get_gamma,matrix(miBeta$coefficients[gamInd],nrow=nbCovs+1),covs=model.matrix(newformula,mhdata),nbStates=nbStates,i=i,j=j)
        se[i,j]<-suppressWarnings(sqrt(dN%*%miBeta$variance[gamInd,gamInd]%*%dN))
        lower[i,j]<-1/(1+exp(-(log(est[i,j]/(1-est[i,j]))-quantSup*(1/(est[i,j]-est[i,j]^2))*se[i,j])))#est[i,j]-quantSup*se[i,j]
        upper[i,j]<-1/(1+exp(-(log(est[i,j]/(1-est[i,j]))+quantSup*(1/(est[i,j]-est[i,j]^2))*se[i,j])))#est[i,j]+quantSup*se[i,j]
      }
    }
    Par$real$gamma <- list(est=est,se=se,lower=lower,upper=upper)
    dimnames(Par$real$gamma$est) <- dimnames(Par$real$gamma$se) <- dimnames(Par$real$gamma$lower) <- dimnames(Par$real$gamma$upper) <- list(m$stateNames,m$stateNames)
  }
  
  # pooled delta estimates
  if(!m$conditions$stationary & nbStates>1){
    nbCovsDelta <- ncol(m$covsDelta)-1
    deltInd <- (length(miBeta$coefficients)-(nbCovsDelta+1)*(nbStates-1)+1):length(miBeta$coefficients)
    delta <- matrix(miBeta$coefficients[deltInd],nrow=nbCovsDelta+1)
    est<-lower<-upper<-se<-matrix(NA,nrow=nrow(m$covsDelta),ncol=nbStates)
    for(j in 1:nrow(m$covsDelta)){
      est[j,] <- get_delta(delta,m$covsDelta[j,,drop=FALSE],1:nbStates)
      for(i in 1:nbStates){
        dN<-numDeriv::grad(get_delta,delta,covsDelta=m$covsDelta[j,,drop=FALSE],i=i)
        se[j,i]<-suppressWarnings(sqrt(dN%*%miBeta$variance[deltInd,deltInd]%*%dN))
        lower[j,i] <- probCI(est[j,i],se[j,i],quantSup,bound="lower")
        upper[j,i] <- probCI(est[j,i],se[j,i],quantSup,bound="upper")
      }
    }
  } else {
    if(nbStates>1){
      covs<-model.matrix(newformula,tempCovs)
      statFun<-function(beta,nbStates,covs,i){
        gamma <- trMatrix_rcpp(nbStates,beta,covs)[,,1]
        tryCatch(solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))[i],error = function(e) {
          "A problem occurred in the calculation of the stationary distribution."})
      }
      delta <- statFun(matrix(miBeta$coefficients[gamInd],nrow=nbCovs+1),nbStates,covs,1:nbStates)
      est <- matrix(delta,nrow=nrow(m$covsDelta),ncol=nbStates,byrow=TRUE)
      lower<-upper<-se<-matrix(NA,nrow=nrow(m$covsDelta),ncol=nbStates)
      for(k in 1:nbStates){
        dN<-numDeriv::grad(statFun,matrix(miBeta$coefficients[gamInd],nrow=nbCovs+1),nbStates=nbStates,covs=covs,i=k)
        se[,k]<-suppressWarnings(sqrt(dN%*%miBeta$variance[gamInd,gamInd]%*%dN))
        lower[,k] <- probCI(est[,k],se[,k],quantSup,bound="lower")
        upper[,k] <- probCI(est[,k],se[,k],quantSup,bound="upper")
      }
    } else {
      est <- matrix(1,nrow(m$covsDelta))
      lower <- upper <- se <- matrix(NA,nrow(m$covsDelta))
    }
  }
  Par$real$delta <- list(est=est,se=se,lower=lower,upper=upper)
  colnames(Par$real$delta$est) <- colnames(Par$real$delta$se) <- colnames(Par$real$delta$lower) <- colnames(Par$real$delta$upper) <- m$stateNames
  rownames(Par$real$delta$est) <- rownames(Par$real$delta$se) <- rownames(Par$real$delta$lower) <- rownames(Par$real$delta$upper) <- paste0("ID:",unique(m$data$ID))
  
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
  }

  if(nbStates>1) {
    Par$timeInStates <- list(est=xbar$timeInStates,se=MI_se$timeInStates,lower=lower$timeInStates,upper=upper$timeInStates)
    names(Par$timeInStates$est) <- m$stateNames
    names(Par$timeInStates$se) <- m$stateNames
    names(Par$timeInStates$lower) <- m$stateNames
    names(Par$timeInStates$upper) <- m$stateNames
    
    Par$states <- states
    
    Par$stateProbs <- list(est=xbar$stateProbs,se=MI_se$stateProbs,lower=lower$stateProbs,upper=upper$stateProbs)
    rownames(Par$stateProbs$est) <- data$ID
    rownames(Par$stateProbs$se) <- data$ID
    rownames(Par$stateProbs$lower) <- data$ID
    rownames(Par$stateProbs$upper) <- data$ID
    colnames(Par$stateProbs$est) <- m$stateNames
    colnames(Par$stateProbs$se) <- m$stateNames
    colnames(Par$stateProbs$lower) <- m$stateNames
    colnames(Par$stateProbs$upper) <- m$stateNames
  }
  
  mh <- im[[1]]
  attr(mh,"class") <- NULL
  mh$mle <- NULL
  mh$mod <- NULL
  mh$CIreal <- NULL
  mh$CIbeta <- NULL
  if(any(ident)) mh$conditions$fullDM <- fullDM
  
  mh$data<-mhdata
  mh$rawCovs<-mhrawCovs

  errorEllipse<-NULL
  if(all(c("x","y") %in% names(mh$data))){
    checkerrs <- lapply(im,function(x) x$data[match(c("x","y"),names(x$data))])
    ident <- !unlist(lapply(checkerrs,function(x) isTRUE(all.equal(x,checkerrs[[1]]))))
    if(any(ident)){
      # calculate location alpha% error ellipses
      cat("Calculating location",paste0(alpha*100,"%"),"error ellipses... ")
      registerDoParallel(cores=ncores)
      errorEllipse<-foreach(i = 1:nrow(mh$data)) %dopar% {
        tmp <- cbind(unlist(lapply(im,function(x) x$data$x[i])),unlist(lapply(im,function(x) x$data$y[i])))
        if(length(unique(tmp[,1]))>1 | length(unique(tmp[,2]))>1)
          ellip <- car::dataEllipse(tmp,levels=alpha,draw=FALSE,segments=100)
        else ellip <- matrix(tmp[1,],101,2,byrow=TRUE)
      }
      stopImplicitCluster()
      cat("DONE\n")
    }
  }
  mh$errorEllipse <- errorEllipse
  mh$Par <- Par
  mh$MIcombine <- miBeta
  
  return(miSum(mh))
}

mi_parm_list<-function(est,se,lower,upper,m){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  dimnames(Par$est) <- dimnames(Par$se) <- dimnames(Par$lower) <- dimnames(Par$upper) <- list(rownames(m),colnames(m))
  Par
}

probCI<-function(x,se,z,bound="lower"){
  if(bound=="lower")
    ci<-boot::inv.logit(boot::logit(x)-z*(1/(x-x^2))*se) 
  else if(bound=="upper")
    ci<-boot::inv.logit(boot::logit(x)+z*(1/(x-x^2))*se) 
  ci
}
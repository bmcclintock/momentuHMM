#' Get starting values for new model from existing \code{momentuHMM} model fit
#' 
#' For nested models, this function will extract starting parameter values (i.e., \code{Par0} in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}) from an existing \code{\link{momentuHMM}} model fit based on the provided arguments for the new model. Any parameters that are not in common between \code{model} and the new model (as specified by the arguments) are set to \code{0}. This function is intended to help users incrementally build and fit more complicated models from simpler nested models (and vice versa).
#' 
#' @param model A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object (as returned by \code{\link{fitHMM}}, \code{\link{MIfitHMM}}, or \code{\link{MIpool}})
#' @param nbStates Number of states in the new model. Default: \code{nbStates=length(model$stateNames)}
#' @param estAngleMean Named list indicating whether or not the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy') are to be estimated in the new model. Default: \code{estAngleMean=model$conditions$estAngleMean}
#' @param circularAngleMean Named list indicating whether circular-linear or circular-circular 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles are to be used in the new model.  See \code{\link{fitHMM}}. Default: \code{circularAngleMean=model$conditions$circularAngleMean}
#' @param formula Regression formula for the transition probability covariates of the new model (see \code{\link{fitHMM}}).  Default: \code{formula=model$conditions$formula}.
#' @param formulaDelta Regression formula for the initial distribution covariates of the new model (see \code{\link{fitHMM}}).  Default: \code{formulaDelta=model$conditions$formulaDelta}.
#' @param mixtures Number of mixtures for the state transition probabilities  (see \code{\link{fitHMM}}). Default: \code{formula=model$conditions$mixtures}.
#' @param DM Named list indicating the design matrices to be used for the probability distribution parameters of each data stream in the new model (see \code{\link{fitHMM}}). Only parameters with design matrix column names that match those in model$conditions$fullDM are extracted, so care must be taken in naming columns if any elements of \code{DM}
#' are specified as matrices instead of formulas. Default: \code{DM=model$conditions$DM}.
#' @param betaRef Numeric vector of length \code{nbStates} indicating the reference elements for the t.p.m. multinomial logit link. Default: \code{formula=model$conditions$betaRef}.
#' @param stateNames Character vector of length \code{nbStates} indicating the names and order of the states in the new model. Default: \code{stateNames=model$stateNames[1:nbStates]}.
#'
#' @return 
#' A named list containing starting values suitable for \code{Par0} and \code{beta0} arguments in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}:
#' \item{Par}{A list of vectors of state-dependent probability distribution parameters for 
#' each data stream specified in \code{model$conditions$dist}}
#' \item{beta}{Matrix of regression coefficients for the transition probabilities}
#' \item{delta}{Initial distribution of the HMM. Only returned if \code{stateNames} has the same membership as the state names for \code{model}}.
#' 
#' All other \code{\link{fitHMM}} (or \code{\link{MIfitHMM}}) model specifications (e.g., \code{dist}, \code{userBounds}, \code{workBounds}, etc.) and \code{data} are assumed to be the same 
#' for \code{model} and the new model (as specified by  \code{estAngleMean}, \code{circularAngleMean}, \code{formula}, \code{formulaDelta}, \code{DM}, and \code{stateNames}).
#'
#' @seealso \code{\link{getPar}}, \code{\link{getParDM}}, \code{\link{fitHMM}}, \code{\link{MIfitHMM}}
#' 
#' @examples 
#' # model is a momentuHMM object, automatically loaded with the package
#' model <- example$m
#' data <- model$data
#' dist <- model$conditions$dist
#' nbStates <- length(model$stateNames)
#' estAngleMean <- model$conditions$estAngleMean
#' 
#' newformula <- ~cov1+cov2
#' Par0 <- getPar0(model,formula=newformula)
#' 
#' \dontrun{
#' newModel <- fitHMM(model$data,dist=dist,nbStates=nbStates,
#'                    Par0=Par0$Par,beta0=Par0$beta,
#'                    formula=newformula,
#'                    estAngleMean=estAngleMean)
#' }
#' 
#' newDM1 <- list(step=list(mean=~cov1,sd=~cov1))
#' Par0 <- getPar0(model,DM=newDM1)
#' 
#' \dontrun{
#' newModel1 <- fitHMM(model$data,dist=dist,nbStates=nbStates,
#'                    Par0=Par0$Par,beta0=Par0$beta,
#'                    formula=model$conditions$formula,
#'                    estAngleMean=estAngleMean,
#'                    DM=newDM1)
#' }
#' 
#' # same model but specify DM for step using matrices
#' newDM2 <- list(step=matrix(c(1,0,0,0,
#'                            "cov1",0,0,0,
#'                            0,1,0,0,
#'                            0,"cov1",0,0,
#'                            0,0,1,0,
#'                            0,0,"cov1",0,
#'                            0,0,0,1,
#'                            0,0,0,"cov1"),nrow=nbStates*2))
#'                            
#' # to be extracted, new design matrix column names must match 
#' # column names of model$conditions$fullDM
#' colnames(newDM2$step)<-paste0(rep(c("mean_","sd_"),each=2*nbStates),
#'                       rep(1:nbStates,each=2),
#'                       rep(c(":(Intercept)",":cov1"),2*nbStates))
#' Par0 <- getPar0(model,DM=newDM2)
#'                       
#' \dontrun{
#' newModel2 <- fitHMM(model$data,dist=dist,nbStates=nbStates,
#'                    Par0=Par0$Par,beta0=Par0$beta,
#'                    formula=model$conditions$formula,
#'                    estAngleMean=estAngleMean,
#'                    DM=newDM2)
#' }
#' 
#' @export
getPar0<-function(model,nbStates=length(model$stateNames),estAngleMean=model$conditions$estAngleMean,circularAngleMean=model$conditions$circularAngleMean,formula=model$conditions$formula,formulaDelta=model$conditions$formulaDelta,mixtures=model$conditions$mixtures,DM=model$conditions$DM,betaRef=model$conditions$betaRef,stateNames=model$stateNames){
  
  if(!is.momentuHMM(model) & !is.miHMM(model) & !is.miSum(model))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(model)) model <- model$miSum
  
  model$conditions$optInd <- numeric() # extra hack needed for bc
  model <- delta_bc(model)
  model$conditions$optInd <- NULL # extra hack needed for bc
  
  if(is.miSum(model)){
    model$mle <- lapply(model$Par$real,function(x) x$est)
    model$mle$beta <- model$Par$beta$beta$est
    model$mle$delta <- model$Par$real$delta$est
    model$mod <- list()
    if(!is.null(model$conditions$recharge)){
      nbRecovs <- ncol(model$g0covs) + ncol(model$reCovs)
      model$mle$g0 <- c(model$Par$beta$g0$est)
      names(model$mle$g0) <- colnames(model$Par$beta$g0$est)
      model$mle$theta <- c(model$Par$beta$theta$est)
      names(model$mle$theta) <- colnames(model$Par$beta$theta$est)
    } else nbRecovs <- 0
    model$mod$estimate <- expandPar(model$MIcombine$coefficients,model$conditions$optInd,unlist(model$conditions$fixPar),model$conditions$wparIndex,model$conditions$betaCons,length(model$stateNames),model$covsDelta,model$conditions$stationary,nrow(model$Par$beta$beta$est)-1,nbRecovs,model$conditions$mixtures)
    model$CIbeta <- model$Par$beta
    model$CIreal <- model$Par$real
  }
  
  if(!is.null(model$mod$hessian) & inherits(model$CIbeta,"error")){
    model$mod$hessian <- NULL
    model$CIbeta <- tryCatch(CIbeta(model),error=function(e) e)
  }
  
  dist<-model$conditions$dist
  Par<-model$mle
  zeroInflation<-model$conditions$zeroInflation
  oneInflation<-model$conditions$oneInflation
  if(is.null(nbStates)) nbStates<-length(model$stateNames)
  if(is.null(estAngleMean)) estAngleMean<-model$conditions$estAngleMean
  if(is.null(circularAngleMean)) circularAngleMean<-model$conditions$circularAngleMean
  if(is.null(formula)) formula<-model$conditions$formula
  if(is.null(formulaDelta)) {
    formDelta <- ~1
  } else formDelta <- formulaDelta
  if(is.null(DM)) DM<-model$conditions$DM
  if(!is.null(stateNames)){
    if(length(stateNames)!=nbStates) stop("stateNames must be of length ",nbStates)
    if(!any(stateNames %in% model$stateNames)) warning("stateNames do not match any model$stateNames")
  }
  else {
    stateNames<-model$stateNames[1:nbStates]
  }
  
  distnames<-names(dist)
  
  if(nbStates<0)
    stop("nbStates should be at least 1.")
  
  if(!is.list(estAngleMean) | is.null(names(estAngleMean))) stop("'estAngleMean' must be a named list")
  for(i in distnames[which(!(dist %in% angledists))]){
    estAngleMean[[i]] <- FALSE
  }
  for(i in distnames){
    if(is.null(estAngleMean[[i]])) estAngleMean[[i]] <- FALSE
    if(!is.logical(estAngleMean[[i]])) stop("estAngleMean$",i," must be logical")
  }
  estAngleMean<-estAngleMean[distnames]
  
  if(!is.list(circularAngleMean) | is.null(names(circularAngleMean))) stop("'circularAngleMean' must be a named list")
  for(i in distnames){
    if(is.null(circularAngleMean[[i]]) | !estAngleMean[[i]]) circularAngleMean[[i]] <- FALSE
    if(!is.logical(circularAngleMean[[i]]) & !is.numeric(circularAngleMean[[i]]) | length(circularAngleMean[[i]])!=1) stop("circularAngleMean$",i," must be logical or numeric")
    if(!isFALSE(circularAngleMean[[i]]) & is.null(DM[[i]])) stop("DM$",i," must be specified when circularAngleMean$",i," = ",circularAngleMean[[i]])
  }
  circularAngleMean<-circularAngleMean[distnames]
  
  if(!is.null(stateNames) & length(stateNames)!=nbStates)
    stop("stateNames must have length ",nbStates)
  
  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,userBounds=NULL)
  parSize <- p$parSize
  
  bounds <- p$bounds
  
  if(is.null(DM)){
    DM <- vector('list',length(distnames))
    names(DM) <- distnames
  } else {
    if(!is.list(DM) | is.null(names(DM))) stop("'DM' must be a named list")
    if(!any(names(DM) %in% distnames)) stop("DM names must include at least one of: ",paste0(distnames,collapse=", "))
  }

  #Par <- list()
  for(i in distnames){
    if(is.null(DM[[i]])) Par[[i]]<-c(t(model$mle[[i]][,which(colnames(model$mle[[i]]) %in% stateNames)]))
  }
  
  DMinputs<-getDM(model$data,DM,dist,nbStates,p$parNames,p$bounds,Par,cons=NULL,workcons=NULL,zeroInflation,oneInflation,circularAngleMean,FALSE)$fullDM
  
  parCount<- lapply(model$conditions$fullDM,ncol)
  for(i in distnames[!unlist(lapply(model$conditions$circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(model$conditions$fullDM[[i]])))))
  }
  parindex <- c(0,cumsum(unlist(parCount)))
  names(parindex) <- c(distnames,"beta")
  
  for(i in distnames){
    #if(!is.null(DM[[i]])) {
      parmNames<-colnames(model$conditions$fullDM[[i]])
      if(!isFALSE(model$conditions$circularAngleMean[[i]])) parmNames <- unique(gsub("cos","",gsub("sin","",parmNames)))
      for(j in 1:length(model$stateNames)){
        for(k in p$parNames[[i]]){
          parmNames<-gsub(paste0(k,"_",j),paste0(k,"_",model$stateNames[j]),parmNames)
        }
      }
      newparmNames<-colnames(DMinputs[[i]])
      if(!isFALSE(circularAngleMean[[i]])) newparmNames <- unique(gsub("cos","",gsub("sin","",newparmNames)))
      tmpPar <- rep(0,length(newparmNames))
      for(j in 1:length(stateNames)){
        for(k in p$parNames[[i]]){
          newparmNames<-gsub(paste0(k,"_",j),paste0(k,"_",stateNames[j]),newparmNames)
        }
      }
      if(any(newparmNames %in% parmNames)){
        if(is.miSum(model))
          tmpPar[match(parmNames,newparmNames,nomatch=0)] <- nw2w((model$CIbeta[[i]]$est-model$conditions$workcons[[i]])^(1/model$conditions$cons[[i]]),model$conditions$workBounds[[i]])[parmNames %in% newparmNames]#model$CIbeta[[i]]$est[parmNames %in% newparmNames]
        else
          tmpPar[match(parmNames,newparmNames,nomatch=0)] <- model$mod$estimate[parindex[[i]]+1:parCount[[i]]][parmNames %in% newparmNames]
      }
      if(!isFALSE(circularAngleMean[[i]])) names(tmpPar) <- unique(gsub("cos","",gsub("sin","",colnames(DMinputs[[i]]))))
      else names(tmpPar) <- colnames(DMinputs[[i]])
      Par[[i]] <- tmpPar
    #}
    if(is.null(DM[[i]])){
      for(j in model$stateNames){
        if(dist[[i]] %in% angledists & !estAngleMean[[i]]){
          Par[[i]][grepl(j,newparmNames)]<-model$mle[[i]][-1,j]
        } else {
          Par[[i]][grepl(j,newparmNames)]<-model$mle[[i]][,j]
        }
        names(Par[[i]])<-paste0(rep(p$parNames[[i]],each=nbStates),"_",1:nbStates)
      }
    }
  }
  
  if(nbStates>1){
    newForm <- newFormulas(formula,nbStates)
    formulaStates <- newForm$formulaStates
    formterms <- newForm$formterms
    newformula <- newForm$newformula
    recharge <- newForm$recharge
    
    betaRow <- rownames(model$mle$beta)
    
    if(!is.null(recharge)){
      model$data$recharge <- rep(0,nrow(model$data))
      betaNames <- colnames(model.matrix(as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge")),model$data))
    } else {
      betaNames <- colnames(model.matrix(newformula,model$data))
    }
    betaNames <- paste0(rep(betaNames,mixtures),"_mix",rep(1:mixtures,each=length(betaNames)))
    if(model$conditions$mixtures==1) betaRow <- paste0(rep(betaRow,model$conditions$mixtures),"_mix",rep(1:model$conditions$mixtures,each=length(betaRow)))
                                                        
    tmpPar <- matrix(0,length(betaNames),nbStates*(nbStates-1))
    #betaRef <- model$conditions$betaRef
    #if(length(model$stateNames)==1) tmpPar[1,] <- rep(-1.5,nbStates*(nbStates-1))
    columns <- NULL
    for(i in 1:nbStates){
      for(j in 1:nbStates) {
        if(betaRef[i]<j)
          columns[(i-1)*nbStates+j-i] <- paste(i,"->",j)
        if(j<betaRef[i])
          columns[(i-1)*(nbStates-1)+j] <- paste(i,"->",j)
      }
    }
    colnames(tmpPar) <- columns
    rownames(tmpPar) <- betaNames
    
    if(length(which(model$stateNames %in% stateNames))>1){
      for(i in match(model$stateNames,stateNames,nomatch=0)){#which(model$stateNames %in% stateNames)){
        for(j in match(model$stateNames,stateNames,nomatch=0)){#which(model$stateNames %in% stateNames)) {
          if(betaRef[i]!=j & i>0 & j>0){
            betaInd <- paste(match(stateNames[i],model$stateNames),"->",match(stateNames[j],model$stateNames))
            if(betaInd %in% colnames(model$mle$beta)){
              if(is.miSum(model))
                tmpPar[betaNames %in% betaRow,paste(i,"->",j)]<-nw2w(model$mle$beta,model$conditions$workBounds$beta)[betaRow %in% betaNames,betaInd]
              else 
                tmpPar[betaNames %in% betaRow,paste(i,"->",j)]<-matrix(model$mod$estimate[parindex[["beta"]]+1:length(model$mle$beta)],nrow(model$mle$beta),ncol(model$mle$beta),dimnames = dimnames(model$mle$beta))[betaRow %in% betaNames,betaInd]
            }
          }
        }
      }
    }
    for(state in 1:(nbStates*(nbStates-1))){
      noBeta<-which(match(colnames(model.matrix(newformula,model$data)),colnames(model.matrix(formulaStates[[state]],model$data)),nomatch=0)==0)
      if(length(noBeta)) tmpPar[noBeta,state] <- 0
    }
    Par$beta<-tmpPar
    if(mixtures==1){
      if(!is.null(recharge)){
        model$data$recharge <- rep(0,nrow(model$data))
        betaNames <- colnames(model.matrix(as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge")),model$data))
      } else {
        betaNames <- colnames(model.matrix(newformula,model$data))
      }
      rownames(Par$beta) <- betaNames
    }
    if(setequal(colnames(model$mle$delta),stateNames) & nbStates>1) {
      if(!length(attr(terms.formula(formDelta),"term.labels")) & is.null(model$conditions$formulaDelta) & is.null(formulaDelta) & model$conditions$mixtures==mixtures){
        Par$delta<-model$mle$delta[1,,drop=FALSE][,match(colnames(model$mle$delta),stateNames)]
        if(mixtures>1){
          Par$delta <- model$mle$delta[seq(1,nrow(model$mle$delta),nrow(model$mle$delta)/mixtures),,drop=FALSE][,match(colnames(model$mle$delta),stateNames)]
          rownames(Par$delta) <- paste0("mix",1:mixtures)
        } else Par$delta <- c(Par$delta)
      } else {
        deltaNames <- colnames(model.matrix(formDelta,model$data))
        nbDeltaCovs <- length(deltaNames)-1
        deltaNames <- paste0(rep(deltaNames,mixtures),"_mix",rep(1:mixtures,each=length(deltaNames)))
        if(model$conditions$mixtures==1) rownames(model$CIbeta$delta$est) <- paste0(rep(rownames(model$CIbeta$delta$est),model$conditions$mixtures),"_mix",rep(1:model$conditions$mixtures,each=length(rownames(model$CIbeta$delta$est))))
        if(all(grepl("(Intercept)",deltaNames)) & is.null(formulaDelta)){
          delta <- matrix(0,mixtures,nbStates)
          for(mix in 1:mixtures){
            tmp <- rep(0,nbStates)
            if(any(deltaNames[(mix-1)*(nbDeltaCovs+1)+1:(nbDeltaCovs+1)] %in% rownames(model$CIbeta$delta$est)))
              tmp<-c(0,model$CIbeta$delta$est[which(rownames(model$CIbeta$delta$est) %in% deltaNames[(mix-1)*(nbDeltaCovs+1)+1:(nbDeltaCovs+1)]),])
            deltaXB <- model$covsDelta[,1,drop=FALSE]%*%matrix(tmp,nrow=1)
            expdelta <- exp(deltaXB)
            tmpdelta <- (expdelta/rowSums(expdelta))[1,,drop=FALSE]
            for(i in which(!is.finite(rowSums(tmpdelta)))){
              tmp <- exp(Brobdingnag::as.brob(deltaXB[i,]))
              tmpdelta[i,] <- as.numeric(tmp/Brobdingnag::sum(tmp))
            }
            delta[(mix-1)*(nbDeltaCovs+1)+1:(nbDeltaCovs+1),] <- tmpdelta
          }
          colnames(delta) <- stateNames
        } else {
          delta <- matrix(0,nrow=length(deltaNames),ncol=nbStates-1,dimnames = list(deltaNames,stateNames[-1]))
          if(length(which(model$stateNames[-1] %in% stateNames[-1])>1)){
            for(i in match(model$stateNames[-1],stateNames[-1],nomatch=0)){#which(model$stateNames %in% stateNames)){
              if(is.miSum(model))
                delta[deltaNames %in% rownames(model$CIbeta$delta$est),i] <- nw2w(model$CIbeta$delta$est,model$conditions$workBounds$delta)[rownames(model$CIbeta$delta$est) %in% deltaNames,match(stateNames[-1][i],model$stateNames[-1])]
              else
                delta[deltaNames %in% rownames(model$CIbeta$delta$est),i] <- matrix(model$mod$estimate[parindex[["beta"]]+length(model$mle$beta)+(model$conditions$mixtures-1)+1:length(model$CIbeta$delta$est)],nrow(model$CIbeta$delta$est),ncol(model$CIbeta$delta$est),dimnames=dimnames(model$CIbeta$delta$est))[rownames(model$CIbeta$delta$est) %in% deltaNames,match(stateNames[-1][i],model$stateNames[-1])]
            }
          }
        }
        if(mixtures==1) delta <- c(delta)
        Par$delta <- delta
      }
    } else {
      Par$delta<-NULL
    }
    if(!is.null(recharge)){
      parmNames <- names(model$mle$g0)
      g0Names <- colnames(model.matrix(recharge$g0,model$data))
      g0 <- rep(0,length(g0Names))
      if(any(g0Names %in% parmNames)){
        if(is.miSum(model))
          g0[g0Names %in% parmNames]<-nw2w(model$mle$g0,model$conditions$workBounds$g0)[parmNames %in% g0Names]
        else 
          g0[g0Names %in% parmNames]<-model$mod$estimate[parindex[["beta"]]+length(model$mle$beta)+length(model$CIbeta$delta$est)+1:length(model$mle$g0)][parmNames %in% g0Names]
      }
      parmNames <- names(model$mle$theta)
      thetaNames <- colnames(model.matrix(recharge$theta,model$data))
      theta <- rep(0,length(thetaNames))
      if(any(thetaNames %in% parmNames)){
        if(is.miSum(model))
          theta[thetaNames %in% parmNames]<-nw2w(model$mle$theta,model$conditions$workBounds$theta)[parmNames %in% thetaNames]
        else 
          theta[thetaNames %in% parmNames]<-model$mod$estimate[parindex[["beta"]]+length(model$mle$beta)+length(model$CIbeta$delta$est)+length(model$mle$g0)+1:length(model$mle$theta)][parmNames %in% thetaNames]
      }
      Par$beta <- list(beta=Par$beta,g0=g0,theta=theta)
    }
    if(mixtures>1){
      if(model$conditions$mixtures==mixtures){
        pie <- model$mle$pi
      } else pie <- rep(1/mixtures,mixtures)
      if(!is.list(Par$beta)) Par$beta <- list(beta=Par$beta,pi=pie)
      else Par$beta$pi <- pie
    }
  } else {
    Par$beta<-NULL
    Par$delta<-NULL
  }
  list(Par=Par[distnames],beta=Par$beta,delta=Par$delta)
}
#' Get starting values for new model from existing \code{momentuHMM} or \code{momentuHierHMM} model fit
#' 
#' For nested models, this function will extract starting parameter values (i.e., \code{Par0} in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}) from an existing \code{\link{momentuHMM}} or \code{momentuHierHMM} model fit based on the provided arguments for the new model. 
#' Any parameters that are not in common between \code{model} and the new model (as specified by the arguments) are set to \code{0}. This function is intended to help users incrementally build and fit more complicated models from simpler nested models (and vice versa).
#' 
#' @param model A \code{\link{momentuHMM}}, \code{\link{momentuHierHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object (as returned by \code{\link{fitHMM}}, \code{\link{MIfitHMM}}, or \code{\link{MIpool}})
#' @param ... further arguments passed to or from other methods
#' @export
getPar0 <- function(model, ...) {
  UseMethod("getPar0")
}

#' @rdname getPar0
#' @method getPar0 default
#' @param nbStates Number of states in the new model. Default: \code{nbStates=length(model$stateNames)}
#' @param estAngleMean Named list indicating whether or not the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy') are to be estimated in the new model. Default: \code{estAngleMean=model$conditions$estAngleMean}
#' @param circularAngleMean Named list indicating whether circular-linear or circular-circular 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles are to be used in the new model.  See \code{\link{fitHMM}}. Default: \code{circularAngleMean=model$conditions$circularAngleMean}
#' @param formula Regression formula for the transition probability covariates of the new model (see \code{\link{fitHMM}}).  Default: \code{formula=model$conditions$formula}.
#' @param formulaDelta Regression formula for the initial distribution covariates of the new model (see \code{\link{fitHMM}}).  Default: \code{formulaDelta=model$conditions$formulaDelta}.
#' @param stationary \code{FALSE} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param mixtures Number of mixtures for the state transition probabilities  (see \code{\link{fitHMM}}). Default: \code{formula=model$conditions$mixtures}.
#' @param formulaPi Regression formula for the mixture distribution probabilities (see \code{\link{fitHMM}}). Default: \code{formula=model$conditions$formulaPi}. 
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
#' \item{delta}{Initial distribution of the HMM. Only returned if \code{stateNames} has the same membership as the state names for \code{model}}
#' 
#' @details 
#' All other \code{\link{fitHMM}} (or \code{\link{MIfitHMM}}) model specifications (e.g., \code{dist}, \code{hierDist}, \code{userBounds}, \code{workBounds}, etc.) and \code{data} are assumed to be the same 
#' for \code{model} and the new model (as specified by  \code{nbStates}, \code{hierStates}, \code{estAngleMean}, \code{circularAngleMean}, \code{formula}, \code{hierFormula}, \code{formulaDelta}, \code{hierFormulaDelta}, \code{DM}, etc.).
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
getPar0.default <- function(model,nbStates=length(model$stateNames),estAngleMean=model$conditions$estAngleMean,circularAngleMean=model$conditions$circularAngleMean,formula=model$conditions$formula,formulaDelta=model$conditions$formulaDelta,stationary=model$conditions$stationary,mixtures=model$conditions$mixtures,formulaPi=model$conditions$formulaPi,DM=model$conditions$DM,betaRef=model$conditions$betaRef,stateNames=model$stateNames,...){
  
  if(!is.momentuHMM(model) & !is.miHMM(model) & !is.miSum(model))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  model$conditions$optInd <- numeric() # extra hack needed for bc
  model <- delta_bc(model)
  model$conditions$optInd <- NULL # extra hack needed for bc
  
  if(!is.null(model$mod$hessian) & inherits(model$CIbeta,"error")){
    model$mod$hessian <- NULL
    model$CIbeta <- tryCatch(CIbeta(model),error=function(e) e)
  }
  
  nbAnimals <- length(unique(model$data$ID))
  dist<-model$conditions$dist
  Par<-model$mle
  zeroInflation<-model$conditions$zeroInflation
  oneInflation<-model$conditions$oneInflation
  if(is.null(nbStates)) nbStates<-length(model$stateNames)
  if(is.null(estAngleMean)) estAngleMean<-model$conditions$estAngleMean
  if(is.null(circularAngleMean)) circularAngleMean<-model$conditions$circularAngleMean
  if(is.null(formula)) formula<-model$conditions$formula
  if(is.null(formulaPi)) {
    formPi <- ~1
  } else formPi <- formulaPi
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
  
  DMinputs<-getDM(model$data,DM,dist,nbStates,p$parNames,p$bounds,Par,zeroInflation,oneInflation,circularAngleMean,FALSE)$fullDM
  
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
          tmpPar[match(parmNames,newparmNames,nomatch=0)] <- nw2w(model$CIbeta[[i]]$est,model$conditions$workBounds[[i]])[parmNames %in% newparmNames]#model$CIbeta[[i]]$est[parmNames %in% newparmNames]
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
    reForm <- formatRecharge(nbStates,formula,model$conditions$betaRef,model$data)
    formulaStates <- reForm$formulaStates
    formterms <- reForm$formterms
    newformula <- reForm$newformula
    recharge <- reForm$recharge
    
    betaRow <- rownames(model$mle$beta)
    
    if(!is.null(recharge)){
      model$data <- cbind(model$data,reForm$newdata)
    }
    betaNames <- colnames(stats::model.matrix(newformula,model$data))
    
    if(((length(betaNames)-1) | (ncol(stats::model.matrix(formDelta,model$data))-1)) & !(nbAnimals==nrow(unique(reForm$covs)) & stationary)) stationary <- FALSE
    
    betaNames <- paste0(rep(betaNames,mixtures),"_mix",rep(1:mixtures,each=length(betaNames)))
    
    if(model$conditions$mixtures==1) betaRow <- paste0(rep(betaRow,model$conditions$mixtures),"_mix",rep(1:model$conditions$mixtures,each=length(betaRow)))
                                                        
    tmpPar <- matrix(0,length(betaNames),nbStates*(nbStates-1))
    #betaRef <- model$conditions$betaRef
    #if(length(model$stateNames)==1) tmpPar[1,] <- rep(-1.5,nbStates*(nbStates-1))
    
    if(length(betaRef)!=nbStates) stop("dimension mismatch between nbStates and betaRef")
    
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
          if(i>0 & j>0 && betaRef[i]!=j){
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
      noBeta<-which(match(colnames(stats::model.matrix(newformula,model$data)),colnames(model.matrix(formulaStates[[state]],model$data)),nomatch=0)==0)
      if(length(noBeta)) tmpPar[noBeta,state] <- 0
    }
    Par$beta<-tmpPar
    if(mixtures==1){
      betaNames <- colnames(stats::model.matrix(newformula,model$data))
      rownames(Par$beta) <- betaNames
    }
    if(setequal(colnames(model$mle$delta),stateNames) & nbStates>1 & !stationary) {
      
      deltaNames <- colnames(stats::model.matrix(formDelta,model$data))
      nbDeltaCovs <- length(deltaNames)-1
      deltaNames <- paste0(rep(deltaNames,mixtures),"_mix",rep(1:mixtures,each=length(deltaNames)))
        
      if(!length(attr(stats::terms.formula(formDelta),"term.labels")) & is.null(model$conditions$formulaDelta) & is.null(formulaDelta) & model$conditions$mixtures==mixtures & !model$conditions$stationary){
        Par$delta<-matrix(getPar(model)$delta,nrow=mixtures)[,match(colnames(model$mle$delta),stateNames)]
      } else {
        if(model$conditions$mixtures==1 & !model$conditions$stationary) rownames(model$CIbeta$delta$est) <- paste0(rep(rownames(model$CIbeta$delta$est),model$conditions$mixtures),"_mix",rep(1:model$conditions$mixtures,each=length(rownames(model$CIbeta$delta$est))))
        if(all(grepl("(Intercept)",deltaNames)) & is.null(formulaDelta)){
          delta <- matrix(0,mixtures,nbStates)
          for(mix in 1:mixtures){
            tmp <- rep(0,nbStates)
            if(!model$conditions$stationary){
              if(any(deltaNames[(mix-1)*(nbDeltaCovs+1)+1:(nbDeltaCovs+1)] %in% rownames(model$CIbeta$delta$est)))
                tmp<-c(0,model$CIbeta$delta$est[which(rownames(model$CIbeta$delta$est) %in% deltaNames[(mix-1)*(nbDeltaCovs+1)+1:(nbDeltaCovs+1)]),])
            }
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
          rownames(delta) <- deltaNames
        } else {
          delta <- matrix(0,nrow=length(deltaNames),ncol=nbStates-1,dimnames = list(deltaNames,stateNames[-1]))
          if(!model$conditions$stationary){
            if(length(which(model$stateNames[-1] %in% stateNames[-1])>1)){
              for(i in match(model$stateNames[-1],stateNames[-1],nomatch=0)){#which(model$stateNames %in% stateNames)){
                if(is.miSum(model))
                  delta[deltaNames %in% rownames(model$CIbeta$delta$est),i] <- nw2w(model$CIbeta$delta$est,model$conditions$workBounds$delta)[rownames(model$CIbeta$delta$est) %in% deltaNames,match(stateNames[-1][i],model$stateNames[-1])]
                else
                  delta[deltaNames %in% rownames(model$CIbeta$delta$est),i] <- matrix(model$mod$estimate[parindex[["beta"]]+length(model$mle$beta)+length(model$CIbeta$pi$est)+1:length(model$CIbeta$delta$est)],nrow(model$CIbeta$delta$est),ncol(model$CIbeta$delta$est),dimnames=dimnames(model$CIbeta$delta$est))[rownames(model$CIbeta$delta$est) %in% deltaNames,match(stateNames[-1][i],model$stateNames[-1])]
              }
            }
          }
        }
        if(mixtures==1) {
          deltaNames <- colnames(stats::model.matrix(formDelta,model$data))
          rownames(delta) <- deltaNames
          if(!nbDeltaCovs) delta <- c(delta)
        }
        Par$delta <- delta
      }
    } else {
      Par$delta<-NULL
      if(nbStates>1 & !stationary){
        warning("delta could not be specified; set manually if so desired")
      }
    }
    if(!is.null(recharge)){
      parmNames <- names(model$mle$g0)
      
      forms <- expandRechargeFormulas(recharge)
      
      g0Names <- colnames(stats::model.matrix(forms$g0,model$data))
      g0 <- rep(0,length(g0Names))
      if(any(g0Names %in% parmNames)){
        if(is.miSum(model))
          g0[g0Names %in% parmNames]<-nw2w(model$mle$g0,model$conditions$workBounds$g0)[parmNames %in% g0Names]
        else 
          g0[g0Names %in% parmNames]<-model$mod$estimate[parindex[["beta"]]+length(model$mle$beta)+length(model$CIbeta$delta$est)+1:length(model$mle$g0)][parmNames %in% g0Names]
      }
      names(g0) <- g0Names
      parmNames <- names(model$mle$theta)
      thetaNames <- colnames(stats::model.matrix(forms$theta,model$data))
      theta <- rep(0,length(thetaNames))
      if(any(thetaNames %in% parmNames)){
        if(is.miSum(model))
          theta[thetaNames %in% parmNames]<-nw2w(model$mle$theta,model$conditions$workBounds$theta)[parmNames %in% thetaNames]
        else 
          theta[thetaNames %in% parmNames]<-model$mod$estimate[parindex[["beta"]]+length(model$mle$beta)+length(model$CIbeta$delta$est)+length(model$mle$g0)+1:length(model$mle$theta)][parmNames %in% thetaNames]
      }
      names(theta) <- thetaNames
      Par$beta <- list(beta=Par$beta,g0=g0,theta=theta)
    }
    if(mixtures>1) {
      if(!length(attr(stats::terms.formula(formPi),"term.labels")) & is.null(model$conditions$formulaPi) & is.null(formulaPi) & model$conditions$mixtures==mixtures){
        Par$pi<-model$mle$pi[1,]
      } else {
        piNames <- colnames(stats::model.matrix(formPi,model$data))
        nbPiCovs <- length(piNames)-1
        if(all(grepl("(Intercept)",piNames)) & is.null(formulaPi)){
          tmp <- rep(0,mixtures)
          if(model$conditions$mixtures==mixtures && any(piNames %in% rownames(model$CIbeta$pi$est)))
            tmp<-c(0,model$CIbeta$pi$est[which(rownames(model$CIbeta$pi$est) %in% piNames),])
          piXB <- model$covsPi[,1,drop=FALSE]%*%matrix(tmp,nrow=1)
          exppi <- exp(piXB)
          tmppi <- (exppi/rowSums(exppi))[1,,drop=FALSE]
          for(i in which(!is.finite(rowSums(tmppi)))){
            tmp <- exp(Brobdingnag::as.brob(piXB[i,]))
            tmppi[i,] <- as.numeric(tmp/Brobdingnag::sum(tmp))
          }
          pie <- tmppi
          colnames(pie) <- paste0("mix",1:mixtures)
        } else {
          pie <- matrix(0,nrow=length(piNames),ncol=mixtures-1,dimnames = list(piNames,paste0("mix",2:mixtures)))
          if(mixtures==model$conditions$mixtures && length(which(piNames %in% rownames(model$CIbeta$pi$est))>1)){
            for(i in match(rownames(model$CIbeta$pi$est),piNames,nomatch=0)){#which(model$stateNames %in% stateNames)){
              if(is.miSum(model))
                pie[piNames[i],] <- nw2w(model$CIbeta$pi$est,model$conditions$workBounds$pi)[piNames[i],]
              else
                pie[piNames[i],] <- matrix(model$mod$estimate[parindex[["beta"]]+length(model$mle$beta)+1:length(model$CIbeta$pi$est)],nrow(model$CIbeta$pi$est),ncol(model$CIbeta$pi$est),dimnames=dimnames(model$CIbeta$pi$est))[piNames[i],]
            }
          }
        }
        Par$pi <- pie
      }
      if(!is.list(Par$beta)) Par$beta <- list(beta=Par$beta,pi=Par$pi)
      else Par$beta$pi <- Par$pi
    }
  } else {
    Par$beta<-NULL
    Par$delta<-NULL
  }
  list(Par=Par[distnames],beta=Par$beta,delta=Par$delta)
}

#' @method getPar0 miHMM
#' 
#' @export
getPar0.miHMM <- function(model,...){
  model <- model$miSum
  getPar0.miSum(model,...)
}

#' @method getPar0 miSum
#' 
#' @export
getPar0.miSum <- function(model,...){

    model <- formatmiSum(model)
    model$CIbeta <- model$Par$beta
    model$CIreal <- model$Par$real
    
    getPar0.default(model,...)

}

#' @rdname getPar0
#' @method getPar0 hierarchical
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states (see \code{\link{fitHMM}}). Default: \code{hierStates=model$conditions$hierStates}.
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy (see \code{\link{fitHMM}}). Default: \code{hierFormula=model$conditions$hierFormula}.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' 
#' @details 
#' Note that for hierarchical models, \code{getPar0} will return hierarchical data.tree objects (\code{hierBeta} and \code{hierDelta}) with the default values for \code{fixPar}, \code{betaCons}, and \code{deltaCons};
#' if hierarchical t.p.m. or initial distribution parameters are subject to constraints, then these must be set manually on the list object returned by \code{getPar0}.
#' 
#' @export
getPar0.hierarchical<-function(model,hierStates=model$conditions$hierStates,estAngleMean=model$conditions$estAngleMean,circularAngleMean=model$conditions$circularAngleMean,hierFormula=model$conditions$hierFormula,hierFormulaDelta=model$conditions$hierFormulaDelta,mixtures=model$conditions$mixtures,formulaPi=model$conditions$formulaPi,DM=model$conditions$DM,...){
  
  inputHierHMM <- formatHierHMM(data=model$data,hierStates=hierStates,hierDist=model$conditions$hierDist,hierFormula=hierFormula,hierFormulaDelta=hierFormulaDelta,mixtures=mixtures)
  nbStates <- inputHierHMM$nbStates
  formula <- inputHierHMM$formula
  formulaDelta <- inputHierHMM$formulaDelta
  betaCons <- inputHierHMM$hBetaCons
  betaRef <- inputHierHMM$betaRef
  deltaCons <- inputHierHMM$hDeltaCons
  fixPar <- inputHierHMM$hFixPar
  stateNames <- inputHierHMM$stateNames
  
  Par0 <- getPar0.default(model,nbStates,estAngleMean,circularAngleMean,formula,formulaDelta,FALSE,mixtures,formulaPi,DM,betaRef,stateNames)
  
  if(is.list(Par0$beta)){
    Pi <- Par0$beta$pi
  } else {
    Pi <- NULL
  }
  hier <- mapHier(Par0$beta,Pi,Par0$delta,hierBeta=NULL,hierDelta=NULL,fixPar,betaCons,deltaCons,hierStates,inputHierHMM$newformula,formulaDelta,inputHierHMM$data,mixtures,recharge=inputHierHMM$recharge,fill=TRUE)
  if(!is.null(inputHierHMM$recharge)){
    for(j in names(inputHierHMM$recharge)){
      hier$hierBeta[[j]]$g0 <- Par0$beta$g0[names(hier$hierBeta[[j]]$g0)]
      hier$hierBeta[[j]]$theta <- Par0$beta$theta[names(hier$hierBeta[[j]]$theta)]
    }
  }
  Par0$beta <- Par0$delta <- NULL
  Par0$hierBeta <- hier$hierBeta
  Par0$hierDelta <- hier$hierDelta
  
  Par0
}

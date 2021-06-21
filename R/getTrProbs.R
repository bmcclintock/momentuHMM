#' Transition probability matrix
#' 
#' Computation of the transition probability matrix for each time step as a function of the covariates and the regression
#' parameters. 
#' 
#' @param data \code{\link{momentuHMM}} object, \code{\link{momentuHierHMM}} object, \code{\link{miSum}} object, \code{\link{miHMM}} object, \code{\link{momentuHMMData}} object, \code{\link{momentuHierHMMData}} object, or a data frame containing the covariate values. 
#' 
#' If a data frame is provided, then either \code{nbStates} must be specified (for a regular HMM) or \code{hierStates} and \code{hierDist}
#' must be specified (for a hierarchical HMM).
#' 
#' @param ... further arguments passed to or from other methods; most are ignored if \code{data} is a \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} object
#' @export
getTrProbs <- function(data, ...){
 UseMethod("getTrProbs")
}

#' @rdname getTrProbs
#' @method getTrProbs default
#' @param nbStates Number of states. Ignored unless \code{data} is a data frame.
#' @param beta Matrix of regression coefficients for the transition probabilities
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the transition probability parameters ('beta' and, for recharge models, 'g0' and 'theta'). \code{workBounds$beta} must be a k x 2 matrix, where k=\code{length(beta)}.
#' The first column pertains to the lower bound and the second column the upper bound. Ignored unless \code{data} is a data frame.
#' @param formula Regression formula for the transition probability covariates. Ignored unless \code{data} is a data frame.
#' @param mixtures Number of mixtures for the state transition probabilities.  Ignored unless \code{data} is a data frame.
#' @param betaRef Indices of reference elements for t.p.m. multinomial logit link. Ignored unless \code{data} is a data frame.
#' @param stateNames Optional character vector of length nbStates indicating state names. Ignored unless \code{data} is a data frame.
#' @param getCI Logical indicating whether to calculate standard errors and logit-transformed confidence intervals based on fitted \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} object. Default: FALSE.
#' @param covIndex Integer vector indicating specific rows of the data to be used in the calculations. This can be useful for reducing unnecessarily long computation times (paricularly when \code{getCI=TRUE}), e.g., when \code{formula} includes factor covariates (such as \code{ID}) but no temporal covariates. Ignored if \code{data} is not a \code{\link{momentuHMM}}, \code{\link{momentuHierHMM}}, \code{\link{miSum}}, or \code{\link{miHMM}} object.
#' @param alpha Significance level of the confidence intervals (if \code{getCI=TRUE}). Default: 0.95 (i.e. 95\% CIs).
#' 
#' @return If \code{mixtures=1}, an array of dimension \code{nbStates} x \code{nbStates} x \code{nrow(data)} containing the t.p.m for each observation in \code{data}.
#' If \code{mixtures>1}, a list of length \code{mixtures}, where each element is an array of dimension \code{nbStates} x \code{nbStates} x \code{nrow(data)} containing the t.p.m for each observation in \code{data}.
#' 
#' If \code{getCI=TRUE} then a list of arrays is returned (with elements \code{est}, \code{se}, \code{lower}, and \code{upper}).
#' 
#' @examples
#' m <- example$m
#' trProbs <- getTrProbs(m)
#' 
#' # equivalent
#' trProbs <- getTrProbs(m$data,nbStates=2,beta=m$mle$beta,formula=m$conditions$formula)
#' 
#' \dontrun{
#' # calculate SEs and 95% CIs
#' trProbsSE <- getTrProbs(m, getCI=TRUE)
#' 
#' # plot estimates and CIs for each state transition
#' par(mfrow=c(2,2))
#' for(i in 1:2){
#'   for(j in 1:2){
#'     plot(trProbsSE$est[i,j,],type="l", 
#'          ylim=c(0,1), ylab=paste(i,"->",j))
#'     arrows(1:dim(trProbsSE$est)[3],
#'            trProbsSE$lower[i,j,],
#'            1:dim(trProbsSE$est)[3],
#'            trProbsSE$upper[i,j,],
#'            length=0.025, angle=90, code=3, col=gray(.5), lwd=1.3)
#'   }
#' }
#' 
#' # limit calculations to first 10 observations
#' trProbsSE_10 <- getTrProbs(m, getCI=TRUE, covIndex=1:10)
#' }
#' 
#' @export
getTrProbs.default <- function(data,nbStates,beta,workBounds=NULL,formula=~1,mixtures=1,betaRef=NULL,stateNames=NULL, getCI=FALSE, covIndex=NULL, alpha = 0.95, ...)
{  
  
  if(!is.momentuHMM(data) & !is.miSum(data) & !is.miHMM(data)){
    hierArgs <- list(...)
    argNames <- names(hierArgs)[which(names(hierArgs) %in% c("hierStates","hierDist","hierBeta","hierFormula"))]
    
    ## check that the data is a momentuHMMData object or valid data frame
    if(!is.momentuHMMData(data)){ 
      if(missing(nbStates)){
        if(all(c("hierStates","hierDist") %in% argNames)){
          return(getTrProbs.hierarchical(data=data,workBounds=workBounds,mixtures=mixtures,covIndex=covIndex,alpha=alpha,...))
        }
      }
      if(!is.data.frame(data)) stop('data must be a data.frame')
    }
    if(!missing(nbStates)){
      if(any(c("hierStates","hierDist","hierFormula","hierBeta") %in% argNames))
        stop("Either nbStates must be specified (for a regular HMM) or hierStates and hierDist must be specified (for a hierarchical HMM)")
    }
    
    if(!is.null(betaRef)){
      if(length(betaRef)!=nbStates) stop("betaRef must be a vector of length ",nbStates)
      if(!is.numeric(betaRef)) stop("betaRef must be a numeric vector")
      if(min(betaRef)<1 | max(betaRef)>nbStates) stop("betaRef elements must be between 1 and ",nbStates)
    } else {
      betaRef <- 1:nbStates
    }
    betaRef <- as.integer(betaRef)
    
    if(is.null(stateNames)){
      for(i in 1:nbStates)
        stateNames[i] <- paste("state",i)
    } else if(length(stateNames)!=nbStates){
      stop("stateNames must have length ",nbStates)
    }
    
    if(!is.formula(formula))
      stop("Check the argument 'formula'.")
    
    if(is.list(beta)){
      if(!is.null(beta$g0)) g0 <- w2wn(beta$g0,workBounds$g0)
      if(!is.null(beta$theta)) theta <- w2wn(beta$theta,workBounds$theta)
      beta <- beta$beta
    }
    
    if(!is.null(covIndex)) {
      if(!is.numeric(covIndex) || any(covIndex<1 | covIndex>nrow(data$data))) stop("covIndex can only include integers between 1 and ",nrow(data))
      data <- data[covIndex,,drop=FALSE]
    }
    
    reForm <- formatRecharge(nbStates,formula,betaRef,data,par=list(g0=g0,theta=theta))
    covs <- reForm$covs
    
    if(!is.matrix(beta)) stop("'beta' must be a matrix")
    if(nrow(beta)!=(ncol(covs)*mixtures) | ncol(beta)!=(nbStates*(nbStates-1))) stop('beta must be a matrix with',ncol(covs)*mixtures,"rows and",nbStates*(nbStates-1),"columns")
      
    if(!is.null(workBounds$beta)){
      if(!is.matrix(workBounds$beta) || (nrow(workBounds$beta)!=length(beta) | ncol(workBounds$beta)!=2)) stop('workBounds$beta must be a matrix with',length(beta),"rows and 2 columns")
    }
    
    getCI <- FALSE
    
  } else {
    
    if(is.miHMM(data)) data <- data$miSum
    
    data <- delta_bc(data)
    
    if(is.miSum(data)){
      data <- formatmiSum(data)
      if(length(data$conditions$optInd)){
        Sigma <- matrix(0,length(data$mod$estimate),length(data$mod$estimate))
        Sigma[(1:length(data$mod$estimate))[-data$conditions$optInd],(1:length(data$mod$estimate))[-data$conditions$optInd]] <- data$MIcombine$variance
      } else {
        Sigma <- data$MIcombine$variance
      }
      data$mod$Sigma <- Sigma
    } else if(is.null(data$mod$hessian) | inherits(data$mod$Sigma,"error")) getCI <- FALSE
    
    if(!is.null(covIndex)) {
      if(!is.numeric(covIndex) || any(covIndex<1 | covIndex>nrow(data$data))) stop("covIndex can only include integers between 1 and ",nrow(data$data))
      data$data <- data$data[covIndex,,drop=FALSE]
    }
    stateNames <- data$stateNames
    nbStates <- length(stateNames)
    beta <- data$mle$beta
    g0 <- data$mle$g0
    theta <- data$mle$theta
    workBounds <- NULL
    formula <- data$conditions$formula
    mixtures <- data$conditions$mixtures
    betaRef <- data$conditions$betaRef
    reForm <- formatRecharge(nbStates,formula,betaRef,data$data,par=data$mle)
    covs <- reForm$covs
    hierRecharge <- reForm$hierRecharge
    nbG0covs <- reForm$nbG0covs
    nbRecovs <- reForm$nbRecovs
    newformula <- reForm$newformula
  }
  chkDots(...)
  wnbeta <- w2wn(beta,workBounds$beta)
  
  trMat <- list()
  
  for(mix in 1:mixtures){
    if(nbStates>1)
      trMat[[mix]] <- trMatrix_rcpp(nbStates,wnbeta[(mix-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],as.matrix(covs),betaRef)
    else
      trMat[[mix]] <- array(1,dim=c(1,1,nrow(data)))
    dimnames(trMat[[mix]]) <- list(stateNames,stateNames,NULL)
    
    if(getCI){
      
      Sigma <- data$mod$Sigma
      
      tmpSplineInputs<-getSplineFormula(reForm$newformula,data$data,cbind(data$data,reForm$newdata))
      desMat <- stats::model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
      
      nbCovs <- ncol(covs) - 1
      gamInd<-(length(data$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)*mixtures+1):(length(data$mod$estimate))-(ncol(data$covsPi)*(mixtures-1))-ifelse(nbRecovs,nbRecovs+1+nbG0covs+1,0)-ncol(data$covsDelta)*(nbStates-1)*(!data$conditions$stationary)*mixtures
      quantSup<-qnorm(1-(1-alpha)/2)
      
      tmpSig <- Sigma[gamInd[unique(c(data$conditions$betaCons))],gamInd[unique(c(data$conditions$betaCons))]]
      se <- lower <- upper <- array(NA,dim=dim(trMat[[mix]]),dimnames=list(stateNames,stateNames,NULL))
      cat("Computing SEs and ",alpha*100,"% CIs",ifelse(mixtures>1,paste0(" for mixture ",mix,"... "),"... "),sep="")
      for(i in 1:nbStates){
        for(j in 1:nbStates){
          for(id in unique(data$data$ID)){
            ind <- which(data$data$ID==id)
            if(!is.null(hierRecharge)){
              tmpSig <- Sigma[c(gamInd[unique(c(data$conditions$betaCons))],length(data$mod$estimate)-(nbRecovs+nbG0covs+1):0),c(gamInd[unique(c(data$conditions$betaCons))],length(data$mod$estimate)-(nbRecovs+nbG0covs+1):0)]
              allCovs <- cbind(data$data,reForm$newdata)
              if(inherits(data,"hierarchical")) class(allCovs) <- append("hierarchical",class(allCovs))
              spl <- split(allCovs[ind,,drop=FALSE],1:nrow(desMat[ind,,drop=FALSE]))
              dN<-t(mapply(function(x) tryCatch(numDeriv::grad(get_TrProbs_recharge,data$mod$estimate[c(gamInd[unique(c(data$conditions$betaCons))],length(data$mod$estimate)-(nbRecovs+nbG0covs+1):0)],covs=spl[[x]],formula=newformula,hierRecharge=hierRecharge,nbStates=nbStates,i=i,j=j,betaRef=data$conditions$betaRef,betaCons=data$conditions$betaCons,workBounds=rbind(data$conditions$workBounds$beta,data$conditions$workBounds$theta),mixture=mix,allCovs=allCovs[ind,,drop=FALSE][1:x,,drop=FALSE]),error=function(e) NA),1:length(spl)))
            } else {
              dN<-t(apply(desMat[ind,,drop=FALSE],1,function(x) tryCatch(numDeriv::grad(get_gamma,data$mod$estimate[gamInd[unique(c(data$conditions$betaCons))]],covs=matrix(x,1,dimnames=list(NULL,names(x))),nbStates=nbStates,i=i,j=j,betaRef=data$conditions$betaRef,betaCons=data$conditions$betaCons,workBounds=data$conditions$workBounds$beta,mixture=mix),error=function(e) NA)))
            }
            se[i,j,ind]<-t(apply(dN,1,function(x) tryCatch(suppressWarnings(sqrt(x%*%tmpSig%*%x)),error=function(e) NA)))
            lower[i,j,ind]<-1/(1+exp(-(log(trMat[[mix]][i,j,ind]/(1-trMat[[mix]][i,j,ind]))-quantSup*(1/(trMat[[mix]][i,j,ind]-trMat[[mix]][i,j,ind]^2))*se[i,j,ind])))#trMat[[mix]][i,j,]-quantSup*se[i,j]
            upper[i,j,ind]<-1/(1+exp(-(log(trMat[[mix]][i,j,ind]/(1-trMat[[mix]][i,j,ind]))+quantSup*(1/(trMat[[mix]][i,j,ind]-trMat[[mix]][i,j,ind]^2))*se[i,j,ind])))#trMat[[mix]][i,j,]+quantSup*se[i,j]
          }
        }
      }
      cat("DONE\n")
      trMat[[mix]] <- list(est=trMat[[mix]],se=se,lower=lower,upper=upper)
    }
  }
  
  if(mixtures==1){
    return(trMat[[1]])
  } else {
    return(trMat)
  }
}

#' @rdname getTrProbs
#' @method getTrProbs hierarchical
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See \code{\link{fitHMM}}.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of regression coefficients for the transition probabilities at each level of the hierarchy, including initial values ('beta'), parameter equality constraints ('betaCons'), fixed parameters ('fixPar'), and working scale bounds ('workBounds'). See details.
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). See \code{\link{fitHMM}}.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). See \code{\link{fitHMM}}.
#' 
#' @return 
#' If a hierarchical HMM structure is provided, then a hierarchical data structure containing the state transition probabilities for each time step at each level of the hierarchy ('gamma') is returned.
#' 
#' @export
getTrProbs.hierarchical <- function(data,hierStates,hierBeta,workBounds=NULL,hierFormula=NULL,mixtures=1,hierDist, getCI=FALSE, covIndex=NULL, alpha = 0.95, ...){
  
  chkDots(...)
  
  installDataTree()
  
  if(is.momentuHierHMM(data) | is.miSum(data) | is.miHMM(data)){
    trProbs <- getTrProbs.default(data,getCI=getCI,covIndex=covIndex,alpha=alpha)
    hierStates <- data$conditions$hierStates
    hierBeta <- data$conditions$hierBeta
    hierFormula <- data$conditions$hierFormula
    mixtures <- data$conditions$mixtures
    hierDist <- data$conditions$hierDist
    if(is.momentuHierHMM(data) && (is.null(data$mod$hessian) | inherits(data$mod$Sigma,"error"))) getCI <- FALSE
    data <- data$data
    if(!is.null(covIndex)) {
      data <- data[covIndex,,drop=FALSE]
    }
  } else {
    inputHierHMM <- formatHierHMM(data,hierStates=hierStates,hierDist=hierDist,hierBeta=hierBeta,hierDelta=NULL,hierFormula=hierFormula,mixtures=mixtures,workBounds=workBounds,checkData=FALSE)
    if(is.list(inputHierHMM$beta)) beta <- inputHierHMM$beta$beta
    else beta <- inputHierHMM$beta
    trProbs <- getTrProbs.default(inputHierHMM$data,inputHierHMM$nbStates,beta,inputHierHMM$workBounds,inputHierHMM$newformula,mixtures,inputHierHMM$betaRef,inputHierHMM$stateNames,covIndex=covIndex)
    getCI <- FALSE
  }
  
  if(mixtures==1) trProbs <- list(trProbs)
  
  beta <- data.tree::Node$new("getTrProbs")
  
  ref1 <- hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==2)
  
  beta$AddChild("level1",gamma=list())

  for(mix in 1:mixtures){
    if(getCI){
      beta$level1$gamma[[mix]] <- list()
      for(i in names(trProbs[[mix]])){
        beta$level1$gamma[[mix]][[i]] <- trProbs[[mix]][[i]][ref1,ref1,which(data$level=="1"),drop=FALSE]
        dimnames(beta$level1$gamma[[mix]][[i]]) <- list(names(ref1),names(ref1),NULL)
      }
    } else {
      beta$level1$gamma[[mix]] <- trProbs[[mix]][ref1,ref1,which(data$level=="1"),drop=FALSE]
      dimnames(beta$level1$gamma[[mix]]) <- list(names(ref1),names(ref1),NULL)
    }
  }
  if(mixtures==1) beta$level1$gamma <- beta$level1$gamma[[1]]
  
  for(j in 2:(hierStates$height-1)){
    
    t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
    names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
    
    beta$AddChild(paste0("level",j))
    
    for(k in names(t)){
      ref <- t[[k]]$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)
      if(!is.null(ref)){
        beta[[paste0("level",j)]]$AddChild(k,gamma=list())
        for(mix in 1:mixtures){
          if(getCI){
            beta[[paste0("level",j)]][[k]]$gamma[[mix]] <- list()
            for(i in names(trProbs[[mix]])){
              beta[[paste0("level",j)]][[k]]$gamma[[mix]][[i]] <- trProbs[[mix]][[i]][ref,ref,which(data$level==j),drop=FALSE]
              dimnames(beta[[paste0("level",j)]][[k]]$gamma[[mix]][[i]]) <- list(names(ref),names(ref),NULL)
            }
          } else {
            beta[[paste0("level",j)]][[k]]$gamma[[mix]] <- trProbs[[mix]][ref,ref,which(data$level==j),drop=FALSE]
            dimnames(beta[[paste0("level",j)]][[k]]$gamma[[mix]]) <- list(names(ref),names(ref),NULL)
          }
        }
        if(mixtures==1) beta[[paste0("level",j)]][[k]]$gamma <- beta[[paste0("level",j)]][[k]]$gamma[[1]]
      }
    }  
  }
  return(beta)
}

get_TrProbs_recharge <- function(beta,covs,formula,hierRecharge,nbStates,i,j,betaRef,betaCons,workBounds=NULL,mixture=1,allCovs){
  
  recharge <- expandRechargeFormulas(hierRecharge)
  
  recovs <- stats::model.matrix(recharge$theta,allCovs)
  g0covs <- stats::model.matrix(recharge$g0,allCovs[1,,drop=FALSE])
  
  tmpBeta <- rep(NA,length(betaCons))
  tmpBeta[unique(c(betaCons))] <- beta[1:(length(beta)-(ncol(recovs)+ncol(g0covs)))]
  tmpBeta <- matrix(tmpBeta[betaCons],nrow(betaCons),ncol(betaCons))
  beta <- w2wn(c(tmpBeta,beta[length(beta)-(ncol(recovs)+ncol(g0covs)-1):0]),workBounds)
  
  g0 <- beta[length(beta)-(ncol(recovs)+ncol(g0covs))+1:ncol(g0covs)] %*% t(g0covs)
  theta <- beta[length(beta)-(ncol(recovs)-1):0]
  
  if(inherits(allCovs,"hierarchical")){
    recLevels <- length(hierRecharge)
    recLevelNames <- names(hierRecharge)
    rechargeNames <- paste0("recharge",gsub("level","",recLevelNames))
    colInd <- lapply(recLevelNames,function(x) which(grepl(paste0("I((level == \"",gsub("level","",x),"\")"),colnames(recovs),fixed=TRUE)))
  } else {
    recLevels <- 1
    rechargeNames <- "recharge"
    colInd <- list(1:ncol(recovs))
  }
  for(iLevel in 1:recLevels){
    covs[,rechargeNames[iLevel]] <- g0 + sum(theta[colInd[[iLevel]]]%*%t(recovs[,colInd[[iLevel]],drop=FALSE])) #covs[,rechargeNames[iLevel],drop=FALSE] + theta[colInd[[iLevel]]]%*%t(recovs[,colInd[[iLevel]],drop=FALSE]) # g0  + theta%*%t(recovs)
  }
  
  newcovs <- stats::model.matrix(formula,covs)
  beta <- matrix(beta[1:(length(beta)-(ncol(recovs)+ncol(g0covs)))],ncol=nbStates*(nbStates-1))
  gamma <- trMatrix_rcpp(nbStates,beta[(mixture-1)*ncol(newcovs)+1:ncol(newcovs),,drop=FALSE],newcovs,betaRef)[,,1]
  gamma[i,j]
}

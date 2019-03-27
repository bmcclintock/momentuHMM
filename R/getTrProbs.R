#' Transition probability matrix
#' 
#' Computation of the transition probability matrix for each time step as a function of the covariates and the regression
#' parameters. 
#' 
#' @param data \code{\link{momentuHMM}} object, \code{\link{momentuHierHMM}} object, \code{\link{momentuHMMData}} object, \code{\link{momentuHierHMMData}} object, or a data frame containing the covariate values. 
#' 
#' If a data frame is provided, then either \code{nbStates} must be specified (for a regular HMM) or \code{hierStates} and \code{hierDist}
#' must be specified (for a hierarchical HMM).
#' 
#' @param ... further arguments passed to or from other methods; ignored if \code{data} is a \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} object
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
#' 
#' @return If \code{mixtures=1}, an array of dimension \code{nbStates} x \code{nbStates} x \code{nrow(data)} containing the t.p.m for each observation in \code{data}.
#' If \code{mixtures>1}, a list of length \code{mixtures}, where each element is an array of dimension \code{nbStates} x \code{nbStates} x \code{nrow(data)} containing the t.p.m for each observation in \code{data}.
#' 
#' @examples
#' m <- example$m
#' trProbs <- getTrProbs(m)
#' 
#' # equivalent
#' trProbs <- getTrProbs(m$data,nbStates=2,beta=m$mle$beta,formula=m$conditions$formula)
#' 
#' @export
getTrProbs.default <- function(data,nbStates,beta,workBounds=NULL,formula=~1,mixtures=1,betaRef=NULL,stateNames=NULL, ...)
{  
  
  if(!is.momentuHMM(data)){
    hierArgs <- list(...)
    argNames <- names(hierArgs)[which(names(hierArgs) %in% c("hierStates","hierDist","hierBeta","hierFormula"))]
    
    ## check that the data is a momentuHMMData object or valid data frame
    if(!is.momentuHMMData(data)){ 
      if(missing(nbStates)){
        if(all(c("hierStates","hierDist") %in% argNames)){
          return(getTrProbs.hierarchical(data=data,hierStates=hierArgs$hierStates,hierBeta=hierArgs$hierBeta,workBounds=workBounds,hierFormula=hierArgs$hierFormula,mixtures=mixtures,hierDist=hierArgs$hierDist))
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
    
    reForm <- formatRecharge(nbStates,formula,data,par=list(g0=g0,theta=theta))
    covs <- reForm$covs
    
    if(!is.matrix(beta)) stop("'beta' must be a matrix")
    if(nrow(beta)!=(ncol(covs)*mixtures) | ncol(beta)!=(nbStates*(nbStates-1))) stop('beta must be a matrix with',ncol(covs)*mixtures,"rows and",nbStates*(nbStates-1),"columns")
      
    if(!is.null(workBounds$beta)){
      if(!is.matrix(workBounds$beta) || (nrow(workBounds$beta)!=length(beta) | ncol(workBounds$beta)!=2)) stop('workBounds$beta must be a matrix with',length(beta),"rows and 2 columns")
    }
  } else {
    stateNames <- data$stateNames
    nbStates <- length(stateNames)
    beta <- data$mle$beta
    g0 <- data$mle$g0
    theta <- data$mle$theta
    workBounds <- NULL
    formula <- data$conditions$formula
    mixtures <- data$conditions$mixtures
    betaRef <- data$conditions$betaRef
    reForm <- formatRecharge(nbStates,formula,data$data,par=list(g0=g0,theta=theta))
    covs <- reForm$covs
  }
  wnbeta <- w2wn(beta,workBounds$beta)
  
  trMat <- list()
  
  for(mix in 1:mixtures){
    if(nbStates>1)
      trMat[[mix]] <- trMatrix_rcpp(nbStates,wnbeta[(mix-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],as.matrix(covs),betaRef)
    else
      trMat[[mix]] <- array(1,dim=c(1,1,nrow(data)))
    dimnames(trMat[[mix]]) <- list(stateNames,stateNames,NULL)
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
getTrProbs.hierarchical <- function(data,hierStates,hierBeta,workBounds=NULL,hierFormula=NULL,mixtures=1,hierDist,...){
  
  if(is.momentuHierHMM(data)){
    trProbs <- getTrProbs.default(data)
    hierStates <- data$conditions$hierStates
    hierBeta <- data$conditions$hierBeta
    hierFormula <- data$conditions$hierFormula
    mixtures <- data$conditions$mixtures
    hierDist <- data$conditions$hierDist
    data <- data$data
  } else {
    inputHierHMM <- formatHierHMM(data,hierStates=hierStates,hierDist=hierDist,hierBeta=hierBeta,hierDelta=NULL,hierFormula=hierFormula,mixtures=mixtures,workBounds=workBounds,checkData=FALSE)
    if(is.list(inputHierHMM$beta)) beta <- inputHierHMM$beta$beta
    else beta <- inputHierHMM$beta
    trProbs <- getTrProbs.default(inputHierHMM$data,inputHierHMM$nbStates,beta,inputHierHMM$workBounds,inputHierHMM$newformula,mixtures,inputHierHMM$betaRef,inputHierHMM$stateNames)
  }
  
  if(mixtures==1) trProbs <- list(trProbs)
  
  beta <- data.tree::Node$new("getTrProbs")
  
  ref1 <- hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2)
  
  beta$AddChild("level1",gamma=list())

  for(mix in 1:mixtures){
    beta$level1$gamma[[mix]] <- trProbs[[mix]][ref1,ref1,which(data$level=="1"),drop=FALSE]
    dimnames(beta$level1$gamma[[mix]]) <- list(names(ref1),names(ref1),NULL)
  }
  if(mixtures==1) beta$level1$gamma <- beta$level1$gamma[[1]]
  
  for(j in 2:(hierStates$height-1)){
    
    t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
    names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
    
    beta$AddChild(paste0("level",j))
    
    for(k in names(t)){
      ref <- t[[k]]$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)
      if(!is.null(ref)){
        beta[[paste0("level",j)]]$AddChild(k,gamma=list())
        for(mix in 1:mixtures){
          beta[[paste0("level",j)]][[k]]$gamma[[mix]] <- trProbs[[mix]][ref,ref,which(data$level==j),drop=FALSE]
          dimnames(beta[[paste0("level",j)]][[k]]$gamma[[mix]]) <- list(names(ref),names(ref),NULL)
        }
        if(mixtures==1) beta[[paste0("level",j)]][[k]]$gamma <- beta[[paste0("level",j)]][[k]]$gamma[[1]]
      }
    }  
  }
  beta
}
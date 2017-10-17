
#' Confidence intervals for working (i.e., beta) parameters
#'
#' Computes the standard errors and confidence intervals on the beta (i.e., working) scale of the data stream probability distribution parameters,
#' as well as for the transition probabilities regression parameters. Working scale depends on the real (i.e., natural) scale of the parameters. For 
#' non-circular distributions or for circular distributions with \code{estAngleMean}=FALSE:
#' 
#' 1) if both lower and upper bounds are finite then logit is the working scale;
#' 2) if lower bound is finite and upper bound is infinite then log is the working scale.
#'
#' For circular distributions with \code{estAngleMean}=TRUE and no constraints imposed by a design matrix (DM) or bounds (userBounds), then the working parameters 
#' are complex functions of both the angle mean and concentrations/sd natural parameters (in this case, it's probably best just to focus on the real parameter
#' estimates!).  However, if constraints are imposed by DM or userBounds on circular distribution parameters with \code{estAngleMean}=TRUE and \code{circularAngleMean}=FALSE:
#' 
#' 1) if the natural bounds are (-pi,pi] then tangent is the working scale, otherwise if both lower and upper bounds are finite then logit is the working scale;
#' 2) if lower bound is finite and upper bound is infinite then log is the working scale.
#' 
#' When circular-circular regression is specified using \code{circularAngleMean}, the working scale 
#' for the mean turning angle is not as easily interpretable, but the 
#' link function is atan2(sin(X)*B,1+cos(X)*B), where X are the angle covariates and B the angle coefficients. 
#' Under this formulation, the reference turning angle is 0 (i.e., movement in the same direction as the previous time step). 
#' In other words, the mean turning angle is zero when the coefficient(s) B=0.
#' 
#' @param m A \code{momentuHMM} object
#' @param alpha Significance level of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#'
#' @return A list of the following objects:
#' \item{...}{List(s) of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the working parameters of the data streams}
#' \item{beta}{List of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the working parameters of the transition probabilities}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' CIbeta(m)
#'
#' @export
#' @importFrom MASS ginv
#' @importFrom utils tail

CIbeta <- function(m,alpha=0.95)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  if(length(m$mod)<=1)
    stop("The given model hasn't been fitted.")

  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1.")

  nbStates <- length(m$stateNames)
  
  dist <- m$conditions$dist
  distnames <- names(dist)
  fullDM <- m$conditions$fullDM
  
  m <- delta_bc(m)

  # identify covariates
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
  covs <- model.matrix(newformula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column


  # inverse of Hessian
  Sigma <- ginv(m$mod$hessian)

  p <- parDef(dist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds)
  bounds <- p$bounds
  #if(!all(unlist(lapply(p$bounds,is.numeric)))){
  #  for(i in distnames){
  #    if(!is.numeric(bounds[[i]])){
  #      bounds[[i]] <- gsub(i,"",bounds[[i]],fixed=TRUE)
  #    }
  #  }
  #}
  
  parindex <- c(0,cumsum(unlist(lapply(fullDM,ncol)))[-length(fullDM)])
  names(parindex) <- distnames
  
  # define appropriate quantile
  quantSup <- qnorm(1-(1-alpha)/2)
  
  wpar <- m$mod$estimate
  
  Par <- list()
  for(i in distnames){
    est <- wpar[parindex[[i]]+1:ncol(fullDM[[i]])]^m$conditions$cons[[i]]+m$conditions$workcons[[i]]
    var <- ((m$conditions$cons[[i]]*(wpar[parindex[[i]]+1:ncol(fullDM[[i]])]^(m$conditions$cons[[i]]-1)))^2)*(diag(Sigma)[parindex[[i]]+1:ncol(fullDM[[i]])])
    # if negative variance, replace by NA
    var[which(var<0)] <- NA
    
    # compute lower and upper for working parameters
    wse <- sqrt(var)
    wlower <- est-quantSup*wse
    wupper <- est+quantSup*wse
    
    Par[[i]]<-beta_parm_list(matrix(est,ncol=length(est),byrow=T),matrix(wse,ncol=length(est),byrow=T),matrix(wlower,ncol=length(est),byrow=T),matrix(wupper,ncol=length(est),byrow=T),m$conditions$fullDM[[i]])

  }

  # group CIs for t.p. coefficients
  if(nbStates>1){
    est <- wpar[tail(cumsum(unlist(lapply(fullDM,ncol))),1)+1:((nbCovs+1)*nbStates*(nbStates-1))]
    var <- diag(Sigma)[tail(cumsum(unlist(lapply(fullDM,ncol))),1)+1:((nbCovs+1)*nbStates*(nbStates-1))]
    
    # if negative variance, replace by NA
    var[which(var<0)] <- NA
    
    wse <- sqrt(var)
    wlower <- est-quantSup*wse
    wupper <- est+quantSup*wse
    Par$beta <- list(est=matrix(est,nrow=1+nbCovs),se=matrix(wse,nrow=1+nbCovs),lower=matrix(wlower,nrow=1+nbCovs),upper=matrix(wupper,nrow=1+nbCovs))
  }

  if(!is.null(m$mle$beta)) {
    rownames(Par$beta$est) <- rownames(m$mle$beta)
    rownames(Par$beta$se) <- rownames(m$mle$beta)
    rownames(Par$beta$lower) <- rownames(m$mle$beta)
    rownames(Par$beta$upper) <- rownames(m$mle$beta)
    colnames(Par$beta$est) <- colnames(m$mle$beta)
    colnames(Par$beta$se) <- colnames(m$mle$beta)
    colnames(Par$beta$lower) <- colnames(m$mle$beta)
    colnames(Par$beta$upper) <- colnames(m$mle$beta)
  }
  
  # group CIs for initial distribution
  if(nbStates>1 & !m$conditions$stationary){
    nbCovsDelta <- ncol(m$covsDelta)-1
    foo <- length(wpar)-(nbCovsDelta+1)*(nbStates-1)+1
    est <- wpar[foo:length(wpar)]
    var <- diag(Sigma)[foo:length(wpar)]
    
    # if negative variance, replace by NA
    var[which(var<0)] <- NA
    
    wse <- sqrt(var)
    wlower <- est-quantSup*wse
    wupper <- est+quantSup*wse
    Par$delta <- list(est=matrix(est,nrow=1+nbCovsDelta,ncol=nbStates-1),se=matrix(wse,nrow=1+nbCovsDelta,ncol=nbStates-1),lower=matrix(wlower,nrow=1+nbCovsDelta,ncol=nbStates-1),upper=matrix(wupper,nrow=1+nbCovsDelta,ncol=nbStates-1))
    rownames(Par$delta$est) <- colnames(m$covsDelta)
    rownames(Par$delta$se) <- colnames(m$covsDelta)
    rownames(Par$delta$lower) <- colnames(m$covsDelta)
    rownames(Par$delta$upper) <- colnames(m$covsDelta)  
    colnames(Par$delta$est) <- m$stateNames[-1]
    colnames(Par$delta$se) <- m$stateNames[-1]
    colnames(Par$delta$lower) <- m$stateNames[-1]
    colnames(Par$delta$upper) <- m$stateNames[-1]  
  } #else if(nbStates==1) {
  #  Par$delta <- list(est=matrix(0,1+nbCovsDelta,ncol=nbStates),se=matrix(NA,1+nbCovsDelta,ncol=nbStates),lower=matrix(NA,1+nbCovsDelta,ncol=nbStates),upper=matrix(NA,1+nbCovsDelta,ncol=nbStates))
  #  rownames(Par$delta$est) <- colnames(m$covsDelta)
  #  rownames(Par$delta$se) <- colnames(m$covsDelta)
  #  rownames(Par$delta$lower) <- colnames(m$covsDelta)
  #  rownames(Par$delta$upper) <- colnames(m$covsDelta)  
  #  colnames(Par$delta$est) <- m$stateNames
  #  colnames(Par$delta$se) <- m$stateNames
  #  colnames(Par$delta$lower) <- m$stateNames
  #  colnames(Par$delta$upper) <- m$stateNames
  #}
  return(Par)
}

beta_parm_list<-function(est,se,lower,upper,m){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  rownames(Par$est) <- rownames(Par$se) <- rownames(Par$lower) <- rownames(Par$upper) <- "[1,]"
  colnames(Par$est) <- colnames(m)
  colnames(Par$se) <- colnames(m)
  colnames(Par$lower) <- colnames(m)
  colnames(Par$upper) <- colnames(m)
  Par
}

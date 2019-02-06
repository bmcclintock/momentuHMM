
#' Check parameter length and order for a \code{\link{fitHierHMM}} (or \code{\link{MIfitHierHMM}}) model
#' 
#' Prints parameters with labels based on \code{DM}, \code{hierFormula}, and/or \code{formulaDelta}.  See \code{\link{fitHierHMM}} for 
#' further argument details.
#' 
#' @param data \code{\link{momentuHierHMMData}} object or a data frame containing the data stream and covariate values
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See details.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). Currently
#' supported distributions are 'bern', 'beta', 'exp', 'gamma', 'lnorm', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. See details.
#' @param Par0 Optional named list containing vectors of state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}.  If \code{Par0} is not provided, then ordered parameter indices are returned. 
#' @param beta0 Optional matrix of regression coefficients for the transition probabilities. If \code{beta0} is not provided, then ordered parameter indices are returned. 
#' @param delta0 Optional values or regression coefficients for the initial distribution of the HMM. If \code{delta0} is not provided, then ordered parameter indices are returned. 
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). 
#' @param circularAngleMean An optional named list indicating whether to use circular-linear or circular-circular
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles. 
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). Default: \code{NULL} (no covariate effects). In addition to allowing standard functions in R formulas
#' (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}), special functions include \code{cosinor(cov,period)} for modeling cyclical patterns, spline functions 
#' (\code{\link[splines]{bs}}, \code{\link[splines]{ns}}, \code{\link[splines2]{bSpline}}, \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}, and \code{\link[splines2]{mSpline}}), 
#' and state- or parameter-specific formulas (see details).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities within a given level of the hierarchy.
#' @param formulaDelta Regression formula for the initial distribution. 
#' @param mixtures Number of mixtures for the state transition probabilities.
#' @param formulaPi Regression formula for the mixture distribution probabilities. 
#' Note that only the covariate values from the first row for each individual ID in \code{data} are used (i.e. time-varying covariates cannot be used for the mixture probabilties).
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. 
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. 
#' @param betaCons Matrix of the same dimension as \code{beta0} composed of integers identifying any equality constraints among the t.p.m. parameters.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. 
#' 
#' @seealso \code{\link{fitHierHMM}}, \code{\link{MIfitHierHMM}}
#' 
#' @export

checkHierPar0 <- function(data,hierStates,hierDist,Par0=NULL,beta0=NULL,delta0=NULL,estAngleMean=NULL,circularAngleMean=NULL,hierFormula=NULL,formulaDelta=~1,mixtures=1,formulaPi=NULL,DM=NULL,userBounds=NULL,workBounds=NULL,betaCons=NULL,fixPar=NULL)
{
  
  ## check that the data is a momentuHierHMMData object or valid data frame
  if(!is.momentuHierHMMData(data)){
    if(!is.data.frame(data)) stop('data must be a data.frame')
    if(is.null(data$ID)) data$ID <- rep(1,nrow(data))
    if(is.null(data$level)) stop("data$level must be specified")
    data <- momentuHierHMMData(data)
  }
  
  inputHierHMM <- formatHierHMM(data,hierStates,hierDist,hierFormula,formulaDelta,mixtures,workBounds,betaCons,fixPar)
  nbStates <- inputHierHMM$nbStates
  dist <- inputHierHMM$dist
  formula <- inputHierHMM$formula
  betaRef <- inputHierHMM$betaRef
  stateNames <- inputHierHMM$stateNames
  
  if(is.null(Par0)){
    
    if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
    distnames<-names(dist)
    if(any(is.na(match(distnames,names(data))))){
      tmpdistnames <- distnames
      for(i in which(is.na(match(distnames,names(data))))){
        if(dist[[distnames[i]]] %in% mvndists){
          if(dist[[distnames[i]]] %in% c("mvnorm2","rw_mvnorm2")){
            tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"))
          } else if(dist[[distnames[i]]] %in% c("mvnorm3","rw_mvnorm3")){
            tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"),paste0(distnames[i],".z"))          
          }
        }
      }
      if(any(is.na(match(tmpdistnames,names(data))))) stop(paste0(tmpdistnames[is.na(match(tmpdistnames,names(data)))],collapse=", ")," not found in data")
    }
    for(i in distnames){
      dist[[i]]<-match.arg(dist[[i]],momentuHMMdists)
      if(!is.null(fixPar[[i]])) stop("fixPar$",i," cannot be specified unless Par0 is specified")
    }
    
    if(is.null(estAngleMean)){
      estAngleMean <- vector('list',length(distnames))
      names(estAngleMean) <- distnames
    } else {
      if(!is.list(estAngleMean) | is.null(names(estAngleMean))) stop("'estAngleMean' must be a named list")
    }
    for(i in distnames[which(!(dist %in% angledists))]){
      estAngleMean[[i]] <- FALSE
    }
    for(i in distnames){
      if(is.null(estAngleMean[[i]])) estAngleMean[[i]] <- FALSE
      if(!is.logical(estAngleMean[[i]])) stop("estAngleMean$",i," must be logical")
    }
    estAngleMean<-estAngleMean[distnames]
    
    if(is.null(circularAngleMean)){
      circularAngleMean <- vector('list',length(distnames))
      names(circularAngleMean) <- distnames
    } else {
      if(!is.list(circularAngleMean) | is.null(names(circularAngleMean))) stop("'circularAngleMean' must be a named list")
    }
    consensus <- vector('list',length(distnames))
    names(consensus) <- distnames
    
    for(i in distnames){
      if(dist[[i]]=="vmConsensus"){
        consensus[[i]] <- TRUE
        estAngleMean[[i]] <- TRUE
        if(is.null(circularAngleMean[[i]]) | isFALSE(circularAngleMean[[i]])) circularAngleMean[[i]] <- TRUE
        if(is.null(DM[[i]])) stop("DM$",i," must be specified when dist$",i,"=vmConsensus")
        #dist[[i]] <- gsub("Consensus","",dist[[i]])
      } else consensus[[i]] <- FALSE
      if(is.null(circularAngleMean[[i]]) | !estAngleMean[[i]]) circularAngleMean[[i]] <- FALSE
      if(!is.logical(circularAngleMean[[i]]) & !is.numeric(circularAngleMean[[i]]) | length(circularAngleMean[[i]])!=1) stop("circularAngleMean$",i," must be logical or numeric")
      if(!isFALSE(circularAngleMean[[i]]) & is.null(DM[[i]])) stop("DM$",i," must be specified when circularAngleMean$",i," = ",circularAngleMean[[i]])
    }
    circularAngleMean<-circularAngleMean[distnames]
    consensus<-consensus[distnames]
    
    # determine whether zero-inflation or one-inflation should be included
    zeroInflation <- oneInflation <- vector('list',length(distnames))
    names(zeroInflation) <- names(oneInflation) <- distnames
    for(i in distnames){
      if(dist[[i]] %in% zeroInflationdists){
        if(length(which(data[[i]]==0))>0) {
          zeroInflation[[i]]<-TRUE
        }
        else 
          zeroInflation[[i]]<-FALSE
      }
      else zeroInflation[[i]]<-FALSE
      if(dist[[i]] %in% oneInflationdists){
        if(length(which(data[[i]]==1))>0) {
          oneInflation[[i]]<-TRUE
        }
        else 
          oneInflation[[i]]<-FALSE
      }
      else oneInflation[[i]]<-FALSE
    }
    p <- parDef(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,userBounds)
    nPar <- lapply(p$bounds,function(x) runif(nrow(x),ifelse(is.finite(x[,1]),x[,1],0),ifelse(is.finite(x[,2]),x[,2],max(x[,1]+1,1.e+6))))
    inputs <- checkInputs(nbStates,dist,nPar,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,cons=NULL,workcons=NULL,stateNames=NULL,checkInflation = TRUE)
    
    # convert RW data
    data <- RWdata(dist,data)
    
    tempCovs <- data[1,]
    DMinputs<-getDM(tempCovs,inputs$DM,inputs$dist,nbStates,inputs$p$parNames,inputs$p$bounds,nPar,inputs$cons,inputs$workcons,zeroInflation,oneInflation,inputs$circularAngleMean,FALSE)
    fullDM <- DMinputs$fullDM
    nc <- meanind <- vector('list',length(distnames))
    names(nc) <- names(meanind) <- distnames
    for(i in distnames){
      nc[[i]] <- apply(fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
      # deal with factors
      if(length(tempCovs)){
        for(j in names(which(unlist(lapply(tempCovs,function(x) inherits(x,"factor")))))){
          if(any(grepl(j,inputs$DM[[i]]))){
            tmpCov <- tempCovs
            for(jj in levels(tempCovs[[j]])){
              tmpCov[[j]] <- factor(jj,levels=levels(tempCovs[[j]]))
              tmpgDM<-getDM(tmpCov,inputs$DM[i],inputs$dist[i],nbStates,inputs$p$parNames[i],inputs$p$bounds[i],nPar[i],inputs$cons[i],inputs$workcons[i],zeroInflation[i],oneInflation[i],inputs$circularAngleMean[i],FALSE)$fullDM[[i]]
              nc[[i]] <- nc[[i]] + apply(tmpgDM,1:2,function(x) !all(unlist(x)==0))
            }
          }
        }
      }
      nc[[i]] <- nc[[i]]>0
      if(!isFALSE(inputs$circularAngleMean[[i]])) {
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
    for(i in distnames){
      bInd <- getboundInd(nc[[i]])
      nPar[[i]] <- nPar[[i]][which(!duplicated(bInd))[bInd]]
    }
    par <- getHierParDM(data,hierStates,hierDist,nPar,zeroInflation,oneInflation,estAngleMean,circularAngleMean,DM,userBounds,workBounds)
    
  } else par <- Par0
  m<-suppressMessages(fitHierHMM(data=data,hierStates=hierStates,hierDist=hierDist,
                             Par0=par,beta0=beta0,delta0=delta0,
                             estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                             hierFormula=hierFormula,formulaDelta=formulaDelta,mixtures=mixtures,formulaPi=formulaPi,
                             DM=DM,userBounds=userBounds,workBounds=workBounds,betaCons=betaCons,fit=FALSE,
                             fixPar=fixPar))
  
  inputs <- checkInputs(nbStates,dist,par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,stateNames,checkInflation = TRUE)
  p <- inputs$p
  DMinputs<-getDM(m$data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,par,inputs$cons,inputs$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,inputs$circularAngleMean)
  
  if(is.null(m$conditions$recharge) & mixtures==1){
    if(is.null(beta0) & nbStates>1) {
      beta <- list(beta=m$mle$beta)
      if(!is.null(fixPar$beta)) stop("fixPar$beta cannot be specified unless beta0 is specified")
    } else {
      beta <- list(beta=beta0)
    }
  } else {
    if(is.null(beta0$beta) & nbStates>1) {
      beta <- list(beta=m$mle$beta)
      if(!is.null(fixPar$beta)) stop("fixPar$beta cannot be specified unless beta0$beta is specified")
    } else {
      beta <- list(beta=beta0$beta)
    }
    if(mixtures>1){
      pie <- beta0$pi
      nbCovsPi <- ncol(m$covsPi)-1
      if(!nbCovsPi){
        if(is.null(beta0$pi)){
          if(!is.null(fixPar$pi)) stop("fixPar$pi cannot be specified unless beta0$pi is specified")
          pie <- matrix(1/mixtures,(nbCovsPi+1),mixtures)
        } else {
          pie <- matrix(pie,(nbCovsPi+1),mixtures)
        }
        pie <- log(pie[-1]/pie[1])
      } else if(is.null(beta0$pi)){
        if(!is.null(fixPar$pi)) stop("fixPar$pi cannot be specified unless beta0$pi is specified")
        pie <- matrix(0,nrow=(nbCovsPi+1),ncol=mixtures-1)
      }
      beta$pi <- pie
    }
    if(!is.null(m$conditions$recharge)){
      if(is.null(beta0$g0) & nbStates>1) {
        beta$g0 <- m$mle$g0
        if(!is.null(fixPar$g0)) stop("fixPar$g0 cannot be specified unless beta0$g0 is specified")
      } else {
        beta$g0 <- beta0$g0
      }
      if(is.null(beta0$theta) & nbStates>1) {
        beta$theta <- m$mle$theta
        if(!is.null(fixPar$theta)) stop("fixPar$theta cannot be specified unless beta0$theta is specified")
      } else {
        beta$theta <- beta0$theta
      }
    }
  }
  
  delta <- delta0
  nbCovsDelta <- ncol(m$covsDelta)-1
  if(!nbCovsDelta){
    if(is.null(delta0)){
      if(!is.null(fixPar$delta)) stop("fixPar$delta cannot be specified unless delta0 is specified")
      delta <- matrix(1/nbStates,(nbCovsDelta+1)*mixtures,nbStates)
    } else {
      delta <- matrix(delta,(nbCovsDelta+1)*mixtures,nbStates)
    }
    delta <- apply(delta,1,function(x) log(x[-1]/x[1]))
  } else if(is.null(delta0)){
    if(!is.null(fixPar$delta)) stop("fixPar$delta cannot be specified unless delta0 is specified")
    delta <- matrix(0,nrow=(nbCovsDelta+1)*mixtures,ncol=nbStates-1)
  }
  
  wpar <- n2w(par,p$bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  
  m$mod <- list()
  m$mod$estimate <- wpar
  m$mod$hessian <- matrix(0,length(wpar),length(wpar))
  
  m$CIreal<-CIreal(m)
  m$CIbeta<-CIbeta(m)
  
  distnames <- names(inputs$dist)
  
  for(i in distnames){
    cat("\n")
    if(is.null(inputs$DM[[i]])) {
      cat(i,                "parameters:\n")      
      cat(rep("-",nchar(i)),"------------\n",sep="")
      tmpPar <- m$mle[[i]]
      if((inputs$dist[[i]] %in% angledists) & !m$conditions$estAngleMean[[i]])
        tmpPar <- tmpPar[-1,,drop=FALSE]
      if(is.null(Par0))
        tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),byrow=TRUE,dimnames=list(rownames(tmpPar),colnames(tmpPar)))
      print(tmpPar)
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat(rep("-",nchar(i)),"----------------------------------\n",sep="")
      tmpPar <- m$CIbeta[[i]]$est
      if(is.null(Par0))
        tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
      print(tmpPar)
    }
  }
  
  if(length(m$stateNames)>1){
    if(!is.null(m$conditions$recharge)){
      cat("\n")
      cat("Initial recharge parameter (g0):\n")
      cat("--------------------------------\n")
      tmpPar <- m$mle$g0
      if(is.null(beta0$g0))
        tmpPar[1:length(tmpPar)] <- 1:length(tmpPar)
      print(tmpPar) 
      cat("\n")
      cat("Recharge function parameters (theta):\n")
      cat("-------------------------------------\n")
      tmpPar <- m$mle$theta
      if(is.null(beta0$theta))
        tmpPar[1:length(tmpPar)] <- 1:length(tmpPar)
      print(tmpPar) 
    }
    
    if(mixtures>1){
      cat("\n")
      if(is.null(m$conditions$formulaPi)) {
        formPi <- ~1
      } else formPi <- m$conditions$formulaPi
      if(!length(attr(terms.formula(formPi),"term.labels")) & is.null(m$conditions$formulaPi)){
        tmpPar <- m$mle$pi[1,]
        rownames(tmpPar)<-NULL
        if(is.null(beta0$pi))
          tmpPar[1:length(tmpPar)] <- 1:length(tmpPar)
        cat("Mixture probabilities (pi):\n")
        cat("---------------------------\n")
        print(tmpPar)
      } else {
        cat("Regression coeffs for the mixture probabilities:\n")
        cat("------------------------------------------------\n")
        tmpPar <- m$CIbeta$pi$est
        if(is.null(beta0$pi))
          tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
        print(tmpPar)
      }
    }
    
    #if(!is.null(m$mle$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities (beta):\n")
    cat("----------------------------------------------------------\n")
    tmpPar <- m$mle$beta
    if(is.null(beta0))
      tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
    else if(!is.null(m$conditions$recharge)){
      if(is.null(beta0$beta)) tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
    }
    if(!is.null(betaCons))
      tmpPar <- matrix(tmpPar[c(betaCons)],nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
    print(tmpPar)
    #}
    
    cat("\n")
    m <- delta_bc(m)
    if(is.null(m$conditions$formulaDelta)) {
      formDelta <- ~1
    } else formDelta <- m$conditions$formulaDelta
    if(!length(attr(terms.formula(formDelta),"term.labels")) & is.null(m$conditions$formulaDelta)){
      tmpPar <- m$mle$delta[1:mixtures,]
      if(mixtures==1) rownames(tmpPar)<-NULL
      else rownames(tmpPar) <- paste0("mix",1:mixtures)
      if(is.null(delta0))
        tmpPar[1:length(tmpPar)] <- 1:length(tmpPar)
      cat("Initial distribution:\n")
      cat("---------------------\n")
      print(tmpPar)
    } else {
      cat("Regression coeffs for the initial distribution:\n")
      cat("---------------------------------------------------\n")
      tmpPar <- m$CIbeta$delta$est
      if(is.null(delta0))
        tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
      print(tmpPar)
    }
    
  }
}
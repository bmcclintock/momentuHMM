
#' Check parameter length and order for a \code{\link{fitHMM}} (or \code{\link{MIfitHMM}}) model
#' 
#' Prints parameters with labels based on \code{DM}, \code{formula}, and/or \code{formulaDelta}.  See \code{\link{fitHMM}} for 
#' further argument details.
#' 
#' @param data \code{\link{momentuHMMData}} object, \code{\link{momentuHierHMMData}} object, or a data frame containing the data stream and covariate values
#' @param ... further arguments passed to or from other methods
#' @export
checkPar0 <- function(data, ...) {
  UseMethod("checkPar0")
}

#' @rdname checkPar0
#' @method checkPar0 default
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. 
#' @param Par0 Optional named list containing vectors of state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}.  If \code{Par0} is not provided, then ordered parameter indices are returned. 
#' @param beta0 Optional matrix of regression coefficients for the transition probabilities. If \code{beta0} is not provided, then ordered parameter indices are returned. 
#' @param delta0 Optional values or regression coefficients for the initial distribution of the HMM. If \code{delta0} is not provided, then ordered parameter indices are returned. 
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). 
#' @param circularAngleMean An optional named list indicating whether to use circular-linear or circular-circular
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles. 
#' @param formula Regression formula for the transition probability covariates. 
#' @param formulaDelta Regression formula for the initial distribution. 
#' @param stationary \code{FALSE} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param mixtures Number of mixtures for the state transition probabilities.
#' @param formulaPi Regression formula for the mixture distribution probabilities. 
#' Note that only the covariate values from the first row for each individual ID in \code{data} are used (i.e. time-varying covariates cannot be used for the mixture probabilties).
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. 
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. 
#' @param betaCons Matrix of the same dimension as \code{beta0} composed of integers identifying any equality constraints among the t.p.m. parameters.
#' @param betaRef Numeric vector of length \code{nbStates} indicating the reference elements for the t.p.m. multinomial logit link.
#' @param deltaCons Matrix of the same dimension as \code{delta0} composed of integers identifying any equality constraints among the initial distribution working scale parameters. Ignored unless a formula is provided in \code{formulaDelta}. 
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. 
#' @param prior A function that returns the log-density of the working scale parameter prior distribution(s). 
#' 
#' @seealso \code{\link{fitHMM}}, \code{\link{MIfitHMM}}
#' 
#' @examples
#' m <- example$m
#' checkPar0(data=m$data, nbStates=2, dist=m$conditions$dist,
#'           estAngleMean = m$conditions$estAngleMean,
#'           formula = m$conditions$formula)
#' 
#' par <- getPar(m)
#' checkPar0(data=m$data, nbStates=2, dist=m$conditions$dist,
#'           estAngleMean = m$conditions$estAngleMean,
#'           formula = m$conditions$formula,
#'           Par0=par$Par, beta0=par$beta, delta0=par$delta)
#'           
#' dummyDat <- data.frame(step=0,angle=0,cov1=0,cov2=0)
#' checkPar0(data=dummyDat, nbStates=2, dist=m$conditions$dist,
#'           estAngleMean = m$conditions$estAngleMean,
#'           formula = m$conditions$formula)
#'
#' \dontrun{
#' simDat <- simData(nbStates=2, dist=m$conditions$dist, Par = par$Par,
#'                   spatialCovs = list(forest=forest),
#'                   centers = matrix(0,1,2),
#'                   nbCovs = 2)
#' checkPar0(data = simDat, nbStates=2, dist=m$conditions$dist,
#'           formula = ~forest,
#'           DM = list(step=list(mean=~cov1, sd=~cov2),
#'                     angle=list(mean=~center1.angle,concentration=~1)),
#'           estAngleMean=list(angle=TRUE),
#'           circularAngleMean=list(angle=TRUE))
#'           
#' par <- list(step=rnorm(8),angle=rnorm(4))
#' beta0 <- matrix(rnorm(4),2,2)
#' delta0 <- c(0.5,0.5)
#' checkPar0(data = simDat, nbStates=2, dist=m$conditions$dist,
#'           Par0 = par, beta0 = beta0, delta0 = delta0,
#'           formula = ~forest,
#'           DM = list(step=list(mean=~cov1, sd=~cov2),
#'                     angle=list(mean=~center1.angle,concentration=~1)),
#'           estAngleMean=list(angle=TRUE),
#'           circularAngleMean=list(angle=TRUE))                
#' }
#' @export
checkPar0.default <- function(data,nbStates,dist,Par0=NULL,beta0=NULL,delta0=NULL,estAngleMean=NULL,circularAngleMean=NULL,formula=~1,formulaDelta=NULL,stationary=FALSE,mixtures=1,formulaPi=NULL,DM=NULL,userBounds=NULL,workBounds=NULL,betaCons=NULL,betaRef=NULL,deltaCons=NULL,stateNames=NULL,fixPar=NULL,prior=NULL,...)
{
  
  hierArgs <- list(...)
  argNames <- names(hierArgs)[which(names(hierArgs) %in% c("hierStates","hierDist","hierBeta","hierDelta","hierFormula","hierFormulaDelta"))]
  
  ## check that the data is a momentuHMMData object or valid data frame
  if(!is.momentuHMMData(data)){
    if(missing(nbStates) & missing(dist)){
      if(all(c("hierStates","hierDist") %in% argNames)){
        if(length(attr(stats::terms.formula(formula),"term.labels"))>0 && is.null(hierArgs$hierFormula)) stop("hierFormula should be specified instead of formula")
        if((!is.null(formulaDelta) && length(attr(stats::terms.formula(formulaDelta),"term.labels"))>0) && is.null(hierArgs$hierFormulaDelta)) stop("hierFormulaDelta should be specified instead of formulaDelta")
        class(data) <- append("hierarchical",class(data))
        return(checkPar0.hierarchical(data,hierStates=hierArgs$hierStates,hierDist=hierArgs$hierDist,Par0,hierBeta=hierArgs$hierBeta,hierDelta=hierArgs$hierDelta,estAngleMean,circularAngleMean,hierFormula=hierArgs$hierFormula,hierFormulaDelta = hierArgs$hierFormulaDelta,mixtures,formulaPi,DM,userBounds,workBounds,fixPar,prior))
      }
    }
    if(!is.data.frame(data)) stop('data must be a data.frame')
    if(is.null(data$ID)) data$ID <- rep(1,nrow(data))
    data <- momentuHMMData(data)
  }
  if(!missing(nbStates) | !missing(dist)){
    if(any(c("hierStates","hierDist") %in% argNames))
      stop("Either nbStates and dist must be specified (for a regular HMM) or hierStates and hierDist must be specified (for a hierarchical HMM)")
  }
  
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
    inputs <- checkInputs(nbStates,dist,nPar,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,stateNames=NULL,checkInflation = TRUE)
    
    # convert RW data
    if(any(unlist(dist) %in% rwdists)){
      data <- RWdata(dist,data,knownStates=NULL)
      knownStates <- data$knownStates
      data$knownStates <- NULL
    }
    
    tempCovs <- data[1,]
    DMinputs<-getDM(tempCovs,inputs$DM,inputs$dist,nbStates,inputs$p$parNames,inputs$p$bounds,nPar,zeroInflation,oneInflation,inputs$circularAngleMean,FALSE)
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
              tmpgDM<-getDM(tmpCov,inputs$DM[i],inputs$dist[i],nbStates,inputs$p$parNames[i],inputs$p$bounds[i],nPar[i],zeroInflation[i],oneInflation[i],inputs$circularAngleMean[i],FALSE)$fullDM[[i]]
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
    par <- getParDM.default(data=data,nbStates=nbStates,dist=dist,Par=nPar,zeroInflation=zeroInflation,oneInflation=oneInflation,estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,DM=DM,userBounds=userBounds,workBounds=workBounds)
    
  } else par <- Par0
  m<-suppressMessages(fitHMM(data=data,nbStates=nbStates,dist=dist,
                             Par0=par,beta0=beta0,delta0=delta0,
                             estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                             formula=formula,formulaDelta=formulaDelta,stationary=stationary,mixtures=mixtures,formulaPi=formulaPi,
                             DM=DM,userBounds=userBounds,workBounds=workBounds,betaCons=betaCons,betaRef=betaRef,deltaCons=deltaCons,fit=FALSE,
                             stateNames=stateNames,fixPar=fixPar,prior=prior))
  
  inputs <- checkInputs(nbStates,dist,par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,DM,m$conditions$userBounds,stateNames,checkInflation = TRUE)
  p <- inputs$p
  par <- par[names(inputs$dist)]
  DMinputs<-getDM(m$data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,par,m$conditions$zeroInflation,m$conditions$oneInflation,inputs$circularAngleMean)
  
  for(i in names(fixPar)){
    if(i %in% names(m$conditions$dist)){
      par[[i]][!is.na(fixPar[[i]])] <- fixPar[[i]][!is.na(fixPar[[i]])]
    }
  }
  
  if(is.null(m$conditions$recharge) & mixtures==1){
    if(is.null(beta0) & nbStates>1) {
      beta <- list(beta=m$mle$beta)
      if(!inherits(data,"hierarchical") && !is.null(fixPar$beta)) stop("fixPar$beta cannot be specified unless beta0 is specified")
    } else {
      beta <- list(beta=beta0)
    }
  } else {
    if(is.null(beta0$beta) & nbStates>1) {
      beta <- list(beta=m$mle$beta)
      if(!inherits(data,"hierarchical") && !is.null(fixPar$beta)) stop("fixPar$beta cannot be specified unless beta0$beta is specified")
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
  if(stationary)
    delta <- NULL
  else if(!nbCovsDelta & !inherits(data,"hierarchical")){
    if(is.null(delta0)){
      if(!inherits(data,"hierarchical") && !is.null(fixPar$delta)) stop("fixPar$delta cannot be specified unless delta0 is specified")
      delta <- matrix(1/nbStates,(nbCovsDelta+1)*mixtures,nbStates)
    } else {
      delta <- matrix(delta,(nbCovsDelta+1)*mixtures,nbStates)
    }
    delta <- apply(delta,1,function(x) log(x[-1]/x[1]))
  } else if(is.null(delta0)){
    if(!inherits(data,"hierarchical") && !is.null(fixPar$delta)) stop("fixPar$delta cannot be specified unless delta0 is specified")
    delta <- matrix(0,nrow=(nbCovsDelta+1)*mixtures,ncol=nbStates-1)
  }
  
  wpar <- n2w(par,p$bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,p$Bndind,inputs$dist)
  
  m$mod$hessian <- matrix(0,length(wpar),length(wpar))
  m$conditions$fit <- TRUE
  
  #m$CIreal<-CIreal(m)
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
    
    if(mixtures>1){
      cat("\n")
      if(is.null(m$conditions$formulaPi)) {
        formPi <- ~1
      } else formPi <- m$conditions$formulaPi
      if(!length(attr(stats::terms.formula(formPi),"term.labels")) & is.null(m$conditions$formulaPi)){
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
    
    if(!inherits(data,"hierarchical")){   
      
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
      
      if(!stationary){
        cat("\n")
        m <- delta_bc(m)
        if(is.null(m$conditions$formulaDelta)) {
          formDelta <- ~1
        } else formDelta <- m$conditions$formulaDelta
        if(!length(attr(stats::terms.formula(formDelta),"term.labels")) & is.null(m$conditions$formulaDelta)){
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
      
    } else return(m)
  }
}

#' @rdname checkPar0
#' @method checkPar0 hierarchical
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See \code{\link{fitHMM}}.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). See \code{\link{fitHMM}}.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the initial matrix of regression coefficients for the transition probabilities at each level of the hierarchy ('beta'). See \code{\link{fitHMM}}.
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the initial values for the initial distribution at each level of the hierarchy ('delta'). See \code{\link{fitHMM}}.
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). See \code{\link{fitHMM}}.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). See \code{\link{fitHMM}}. Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' 
#' @export
checkPar0.hierarchical <- function(data,hierStates,hierDist,Par0=NULL,hierBeta=NULL,hierDelta=NULL,estAngleMean=NULL,circularAngleMean=NULL,hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,formulaPi=NULL,DM=NULL,userBounds=NULL,workBounds=NULL,betaCons=NULL,deltaCons=NULL,fixPar=NULL,prior=NULL,...)
{
  
  ## check that the data is a momentuHierHMMData object or valid data frame
  if(!is.momentuHierHMMData(data)){
    if(!is.data.frame(data)) stop('data must be a data.frame')
    if(is.null(data$ID)) data$ID <- rep(1,nrow(data))
    if(is.null(data$level)) stop("data$level must be specified")
    data <- momentuHierHMMData(data)
  }
  
  if(is.null(hierBeta) & !is.null(fixPar$beta)) stop("fixPar$beta cannot be specified unless hierBeta is specified")
  if(is.null(hierDelta) & !is.null(fixPar$delta)) stop("fixPar$delta cannot be specified unless hierDelta is specified")
  
  inputHierHMM <- formatHierHMM(data,hierStates,hierDist,hierBeta,hierDelta,hierFormula,hierFormulaDelta,mixtures,workBounds,betaCons,deltaCons,fixPar)
  nbStates <- inputHierHMM$nbStates
  dist <- inputHierHMM$dist
  beta0 <- inputHierHMM$beta
  delta0 <- inputHierHMM$delta
  formula <- inputHierHMM$formula
  formulaDelta <- inputHierHMM$formulaDelta
  workBounds <- inputHierHMM$workBounds
  betaCons <- inputHierHMM$betaCons
  deltaCons <- inputHierHMM$deltaCons
  betaRef <- inputHierHMM$betaRef
  stateNames <- inputHierHMM$stateNames
  fixPar <- inputHierHMM$fixPar
  recharge <- inputHierHMM$recharge
  
  m <- checkPar0.default(data=data,nbStates=nbStates,dist=dist,Par0=Par0,beta0=beta0,delta0=delta0,estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,formula=formula,formulaDelta=formulaDelta,stationary=FALSE,mixtures=mixtures,formulaPi=formulaPi,DM=DM,userBounds=userBounds,workBounds=workBounds,betaCons=betaCons,betaRef=betaRef,deltaCons=deltaCons,stateNames=stateNames,fixPar=fixPar,prior=prior)
  
  hier <- mapHier(list(beta=m$mle$beta,g0=m$mle$g0,theta=m$mle$theta),m$mle$pi,m$CIbeta$delta$est,hierBeta=hierBeta,hierDelta=hierDelta,inputHierHMM$hFixPar,inputHierHMM$hBetaCons,inputHierHMM$hDeltaCons,hierStates,inputHierHMM$newformula,m$conditions$formulaDelta,inputHierHMM$data,m$conditions$mixtures,recharge,fill=TRUE)
  
  if(!is.null(recharge)){
    g0 <- m$mle$g0
    theta <- m$mle$theta
    for(j in names(recharge)){
      hier$hierBeta[[j]]$g0 <- g0[names(hier$hierBeta[[j]]$g0)]
      hier$hierBeta[[j]]$theta <- theta[names(hier$hierBeta[[j]]$theta)]
    }
  }
  
  if(!is.list(hier$hierBeta)){
    beta0 <- list(beta=hier$hierBeta)
  } else {
    beta0 <- hier$hierBeta
  }
  delta0 <- hier$hierDelta
  
  cat("\n")
  cat("----------------------------------------------------------\n")
  cat("Regression coeffs for the transition probabilities (beta):\n")
  cat("----------------------------------------------------------\n")
  for(j in 1:(hierStates$height-1)){
  cat("------------------------- ",paste0("level",j)," -----------------------\n")
    if(j>1){
      for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
        tmpPar <- beta0$beta[[paste0("level",j)]][[jj]]$beta
        tmpCons <- beta0$beta[[paste0("level",j)]][[jj]]$betaCons
        if(is.null(hierBeta))
          tmpPar <- tmpCons
        print(tmpPar)
        cat("\n")
      }
    } else {
      tmpPar <- beta0$beta[[paste0("level",j)]]$beta
      tmpCons <- beta0$beta[[paste0("level",j)]]$betaCons
      if(is.null(hierBeta))
        tmpPar <- tmpCons
      print(tmpPar)
      cat("\n")
    }
  }
  cat("----------------------------------------------------------\n")
  
  cat("\n")
  cat("----------------------------------------------------------\n")
  cat("Regression coeffs for the initial distribution (delta):\n")
  cat("----------------------------------------------------------\n")
  for(j in 1:(hierStates$height-1)){
    cat("------------------------ ",paste0("level",j)," ------------------------\n")
    if(j>1){
      for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
        tmpPar <- delta0[[paste0("level",j)]][[jj]]$delta
        tmpCons <- delta0[[paste0("level",j)]][[jj]]$deltaCons
        if(is.null(hierDelta))
          tmpPar <- tmpCons
        print(tmpPar)
        cat("\n")
      }
    } else {
      tmpPar <- delta0[[paste0("level",j)]]$delta
      tmpCons <- delta0[[paste0("level",j)]]$deltaCons
      if(is.null(hierDelta))
        tmpPar <- tmpCons
      print(tmpPar)
      cat("\n")
    }
  }
  cat("----------------------------------------------------------\n")
  
  if(!is.null(recharge)){
    cat("\n\n")
    cat("----------------------------------------------------------\n")
    cat("Initial recharge parameter (g0):\n")
    cat("----------------------------------------------------------\n")
    for(j in names(recharge)){
      cat("------------------------- ",j," -----------------------\n")
      tmpPar <- hier$hierBeta[[j]]$g0
      if(is.null(hierBeta))
        tmpPar[1:length(tmpPar)] <- 1:length(tmpPar)
      print(tmpPar) 
    }
    cat("----------------------------------------------------------\n")
    cat("\n")
    cat("----------------------------------------------------------\n")
    cat("Recharge function parameters (theta):\n")
    cat("----------------------------------------------------------\n")
    for(j in names(recharge)){
      cat("------------------------- ",j," -----------------------\n")
      tmpPar <- hier$hierBeta[[j]]$theta
      if(is.null(hierBeta))
        tmpPar[1:length(tmpPar)] <- 1:length(tmpPar)
      print(tmpPar) 
    }
    cat("----------------------------------------------------------\n")
  }

}
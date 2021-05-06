#' Mixture probabilities
#'
#' For a fitted model, this function computes the probability of each individual being in a particular mixture
#'
#' @param m \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} object
#' @param getCI Logical indicating whether to calculate standard errors and logit-transformed confidence intervals for fitted \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} object. Default: FALSE.
#' @param alpha Significance level of the confidence intervals (if \code{getCI=TRUE}). Default: 0.95 (i.e. 95\% CIs).
#' 
#' @return The matrix of individual mixture probabilities, with element [i,j] the probability
#' of individual i being in mixture j
#' 
#' @details
#' When \code{getCI=TRUE}, it can take a while for large data sets and/or a large number of mixtures because the model likelihood for each individual must be repeatedly evaluated in order to numerically approximate the SEs.
#'
#' @examples
#' \dontrun{
#' nObs <- 100
#' nbAnimals <- 20
#' dist <- list(step="gamma",angle="vm")
#' Par <- list(step=c(100,1000,50,100),angle=c(0,0,0.1,2))
#' 
#' # create sex covariate
#' cov <- data.frame(sex=factor(rep(c("F","M"),each=nObs*nbAnimals/2)))
#' formulaPi <- ~ sex + 0
#' 
#' # Females more likely in mixture 1, males more likely in mixture 2
#' beta <- list(beta=matrix(c(-1.5,-0.5,-1.5,-3),2,2),
#'              pi=matrix(c(-2,2),2,1,dimnames=list(c("sexF","sexM"),"mix2"))) 
#' 
#' data.mix<-simData(nbAnimals=nbAnimals,obsPerAnimal=nObs,nbStates=2,dist=dist,Par=Par,
#'                   beta=beta,formulaPi=formulaPi,mixtures=2,covs=cov) 
#' 
#' Par0 <- list(step=Par$step, angle=Par$angle[3:4])   
#' m.mix <- fitHMM(data.mix, nbStates=2, dist=dist, Par0 = Par0, 
#'                 beta0=beta,formulaPi=formulaPi,mixtures=2)
#'                 
#' mixProbs <- mixtureProbs(m.mix, getCI=TRUE)
#' }
#' @references
#' Maruotti, A., and T. Ryden. 2009. A semiparametric approach to hidden Markov models under longitudinal observations. Statistics and Computing 19: 381-393.
#'
#' @export
mixtureProbs <- function(m, getCI=FALSE, alpha = 0.95){
  
  if(!is.momentuHMM(m) & !is.momentuHierHMM(m))# & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM or momentuHierHMM object")
  
  #if(is.miHMM(m)) m <- m$miSum
  
  m <- delta_bc(m)
  
  #if(is.miSum(m)){
  #  m <- formatmiSum(m)
  #  getCI <- FALSE
  #} else 
  if(is.null(m$mod$hessian) | inherits(m$mod$Sigma,"error")) getCI <- FALSE
  
  nbAnimals <- length(unique(m$data$ID))
  mixtures <- m$conditions$mixtures
  
  if(mixtures==1) getCI <- FALSE
    #stop("No mixtures to assign probabilities (mixtures=1)")
  
  #if(mixtures>1) pie <- m$mle[["pi"]]
  #else pie <- matrix(1,nbAnimals,1)
  
  quantSup<-qnorm(1-(1-alpha)/2)
    
  est <- se <- lower <- upper <- matrix(NA,nbAnimals,mixtures)
  colnames(est) <- colnames(se) <- colnames(lower) <- colnames(upper) <- paste0("mix",1:mixtures)
  rownames(est) <- rownames(se) <- rownames(lower) <- rownames(upper) <- paste0("ID:",unique(m$data$ID))
  for(k in 1:nbAnimals){
    #pInd <- which(mapply(function(x) isTRUE(all.equal(x,0)),pie[k,]))
    ind <- which(m$data$ID==unique(m$data$ID)[k])
    tmp <- m
    tmp$data <- m$data[ind,]
    tmp$covsDelta <- m$covsDelta[k,,drop=FALSE]
    tmp$covsPi <- m$covsPi[k,,drop=FALSE]
    tmp$conditions$knownStates <- rep(NA,nrow(tmp$data))
    tmp$conditions$fullDM <- lapply(m$conditions$fullDM,function(y) matrix(mapply(function(x){if(length(x)>1){return(x[ind])} else {return(x)}},y,SIMPLIFY = TRUE),nrow(y),ncol(y),dimnames=dimnames(y)))
    est[k,] <- get_mixProbs(tmp$mod$wpar,mod=tmp,mixture=1:mixtures)
    #est[k,pInd] <- 0
    if(getCI){
      cat("\rComputing SEs and ",alpha*100,"% CIs for individual ",unique(m$data$ID)[k],"... ",sep="")
      for(mix in 1:mixtures){
        dN<-t(tryCatch(numDeriv::grad(get_mixProbs,tmp$mod$wpar,mod=tmp,mixture=mix),error=function(e) NA))
        se[k,mix] <- sqrt(dN %*% m$mod$Sigma %*% t(dN))
        lower[k,mix] <- 1/(1+exp(-(log(est[k,mix]/(1-est[k,mix]))-quantSup*(1/(est[k,mix]-est[k,mix]^2))*se[k,mix])))
        upper[k,mix] <- 1/(1+exp(-(log(est[k,mix]/(1-est[k,mix]))+quantSup*(1/(est[k,mix]-est[k,mix]^2))*se[k,mix])))
      }
    }
  }
  if(getCI) {
    cat("DONE\n")
    out <- list(est=est,se=se,lower=lower,upper=upper)
  } else out <- est
  return(out)
}

get_mixProbs <- function(optPar,mod,mixture){
    
  dist <- mod$conditions$dist
  nbStates <- length(mod$stateNames)
  data <- mod$data
  mixtures <- mod$conditions$mixtures
  formula <- mod$conditions$formula
  knownStates <- mod$conditions$knownStates
  
  newForm <- newFormulas(formula,nbStates,mod$conditions$betaRef,hierarchical = TRUE)
  formulaStates <- newForm$formulaStates
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  # build design matrix for t.p.m.
  covs <- stats::model.matrix(newformula,data)
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  # build design matrix for recharge model
  if(!is.null(recharge)){
    reForm <- formatRecharge(nbStates,formula,mod$conditions$betaRef,data=data)
    formulaStates <- reForm$formulaStates
    newformula <- reForm$newformula
    recharge <- reForm$recharge
    hierRecharge <- reForm$hierRecharge
    newdata <- reForm$newdata
    g0covs <- reForm$g0covs
    nbG0covs <- ncol(g0covs)-1
    recovs <- reForm$recovs
    nbRecovs <- ncol(recovs)-1
    covs <- reForm$covs
    nbCovs <- ncol(covs)-1
    recovsCol <- get_all_vars(recharge$theta,data)#rownames(attr(stats::terms(formula),"factors"))#attr(stats::terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
    if(!all(names(recovsCol) %in% names(data))){
      recovsCol <- recovsCol[,names(recovsCol) %in% names(data),drop=FALSE]
    }
    g0covsCol <- get_all_vars(recharge$g0,data)#rownames(attr(stats::terms(formula),"factors"))#attr(stats::terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
    if(!all(names(g0covsCol) %in% names(data))){
      g0covsCol <- g0covsCol[,names(g0covsCol) %in% names(data),drop=FALSE]
    }
  } else {
    nbRecovs <- 0
    nbG0covs <- 0
    g0covs <- g0covsCol <- NULL
    recovs <- recovsCol <- NULL
    newdata <- NULL
  }
  
  # check arguments
  distnames<-names(dist)
  
  #covs <- stats::model.matrix(formula,data)
  nbCovs <- ncol(covs)-1
  if(!is.null(recovs)) {
    nbRecovs <- ncol(recovs)-1
    nbG0covs <- ncol(g0covs)-1
  } else nbRecovs <- nbG0covs <- 0
  
  ncmean <- get_ncmean(distnames,mod$conditions$fullDM,mod$conditions$circularAngleMean,nbStates)
  nc <- ncmean$nc
  meanind <- ncmean$meanind
  
  consensus <- vector('list',length(distnames))
  names(consensus) <- distnames
  for(i in distnames){
    consensus[[i]] <- (dist[[i]]=="vmConsensus")
  }
  dist <- lapply(dist,function(x) gsub("Consensus","",x))
  dist <- lapply(dist,function(x) ifelse(grepl("cat",x),"cat",x))
  
  wpar <- expandPar(optPar,mod$conditions$optInd,mod$mod$estimate,mod$conditions$wparIndex,mod$conditions$betaCons,mod$conditions$deltaCons,nbStates,ncol(mod$covsDelta)-1,mod$conditions$stationary,nbCovs,nbRecovs+nbG0covs,mixtures,ncol(mod$covsPi)-1)
  par <- w2n(wpar,mod$conditions$bounds,lapply(mod$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,mod$conditions$estAngleMean,mod$conditions$circularAngleMean,consensus,mod$conditions$stationary,mod$conditions$fullDM,mod$conditions$DMind,nrow(mod$data),dist,mod$conditions$Bndind,nc,meanind,mod$covsDelta,mod$conditions$workBounds,mod$covsPi)
  
  if(nbRecovs){
    for(i in 1:length(unique(data$ID))){
      idInd <- which(data$ID==unique(data$ID)[i])
      if(inherits(data,"hierarchical")) {
        recLevels <- length(recharge)
        recLevelNames <- names(recharge)
        rechargeNames <- paste0("recharge",gsub("level","",recLevelNames))
        colInd <- lapply(recLevelNames,function(x) which(grepl(paste0("I((level == \"",gsub("level","",x),"\")"),colnames(recovs),fixed=TRUE)))
      } else {
        recLevels <- 1
        rechargeNames <- "recharge"
        colInd <- list(1:ncol(recovs))
      }
      for(iLevel in 1:recLevels){
        g0 <- par$g0 %*% t(g0covs[(i-1)*recLevels+iLevel,,drop=FALSE])
        theta <- par$theta
        covs[idInd,grepl(rechargeNames[iLevel],colnames(covs))] <- cumsum(c(g0,theta[colInd[[iLevel]]]%*%t(recovs[idInd[-length(idInd)],colInd[[iLevel]]])))
      }
    }
  }
  
  if(is.null(knownStates)) knownStates <- -1
  else knownStates[which(is.na(knownStates))] <- 0
  
  # NULL arguments don't suit C++
  if(any(unlist(lapply(dist,is.null)))){
    par[which(unlist(lapply(dist,is.null)))]<-matrix(NA)
  }
  
  if(mod$conditions$stationary)
    par$delta <- c(NA)
  if(nbStates==1) {
    par$beta <- matrix(NA)
    #par[["pi"]] <- c(NA)
    par$delta <- c(NA)
    par[distnames] <- lapply(par[distnames],as.matrix)
  }
  
  nbCovsDelta <- ncol(mod$covsDelta)-1
    
  beta <- par$beta
  delta <- par$delta
  pie <- par[["pi"]]
  pInd <- which(mapply(function(x) isTRUE(all.equal(x,0)), pie))
  if (length(pInd)) {
    pie[pInd] <- 1e-100
    pie[-pInd] <- pie[-pInd] - (1e-100 * length(pInd))/(ncol(pie) - length(pInd))
  }
  par[["pi"]] <- matrix(1,1,1)
  
  mixProbs <- lnum <- la <- numeric(mixtures)
  for(mix in 1:mixtures){
    par$beta <- beta[(mix-1)*(nbCovs+1)+1:(nbCovs+1),,drop=FALSE]
    if(!mod$conditions$stationary) par$delta <- delta[mix,,drop=FALSE]
    la[mix] <- nLogLike_rcpp(nbStates,as.matrix(covs),data,names(dist),dist,
                             par,
                             1,mod$conditions$zeroInflation,mod$conditions$oneInflation,mod$conditions$stationary,knownStates,mod$conditions$betaRef,1)
    if(!is.null(mod$prior)) la[mix] <- la[mix] - mod$prior(wpar)
    c <- max(-la[mix]+log(pie[mix]))
    lnum[mix] <- c + log(sum(exp(-la[mix]+log(pie[mix])-c)))  
  }
  c <- max(lnum)
  mixProbs <- exp(lnum - c - log(sum(exp(lnum-c))))
  #mixProbs[pInd] <- 0
  return(mixProbs[mixture])
}

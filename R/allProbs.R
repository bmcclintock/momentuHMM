
#' Matrix of all probabilities
#'
#' Used in functions \code{\link{viterbi}}, \code{\link{logAlpha}}, \code{\link{logBeta}}.
#'
#' @param m Object \code{\link{momentuHMM}} or \code{\link{miSum}}.
#'
#' @return Matrix of all probabilities.
#'
#' @examples
#' \dontrun{
#' P <- momentuHMM:::allProbs(m=example$m)
#' }
#' @importFrom LaplacesDemon dbern

allProbs <- function(m)
{
  
  if(!is.momentuHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM) or miSum object (as returned by MIpool)")
  
  m <- delta_bc(m)
  
  data <- m$data
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  zeroInflation <- m$conditions$zeroInflation
  oneInflation <- m$conditions$oneInflation
  nbObs <- nrow(data)
  
  if(is.miSum(m)){

    Par <- lapply(m$Par$real,function(x) x$est)
    for(i in distnames){
      if(!is.null(m$conditions$DM[[i]]))
        Par[[i]] <- m$Par$beta[[i]]$est
      else if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]])
        Par[[i]] <- Par[[i]][-1,]
      
      m$conditions$cons[[i]]<-rep(1,length(m$conditions$cons[[i]]))
      m$conditions$workcons[[i]]<-rep(0,length(m$conditions$workcons[[i]]))
      m$conditions$workBounds[[i]]<-matrix(c(-Inf,Inf),nrow(m$conditions$workBounds[[i]]),2,byrow=TRUE)
    }
    
    Par<-lapply(Par,function(x) c(t(x)))
    Par<-Par[distnames]
    beta <- m$Par$beta$beta$est
    delta <- m$Par$real$delta$est
    if(!is.null(beta)) m$conditions$workBounds$beta<-matrix(c(-Inf,Inf),length(beta),2,byrow=TRUE)
    if(!is.null(m$Par$beta$delta$est)) m$conditions$workBounds$delta<-matrix(c(-Inf,Inf),length(m$Par$beta$delta$est),2,byrow=TRUE)
    
    g0 <- c(m$Par$beta$g0$est)
    theta <- c(m$Par$beta$theta$est)
    if(!is.null(g0)) m$conditions$workBounds$g0<-matrix(c(-Inf,Inf),length(g0),2,byrow=TRUE)
    if(!is.null(theta)) m$conditions$workBounds$theta<-matrix(c(-Inf,Inf),length(theta),2,byrow=TRUE)
    
    inputs <- checkInputs(nbStates,dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
    p <- inputs$p
    DMinputs<-getDM(data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
    m$conditions$fullDM <- DMinputs$fullDM
    m$mod$estimate <- n2w(Par,p$bounds,list(beta=beta,g0=g0,theta=theta),m$Par$beta$delta$est,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
    g0 <- m$mle$g0
    theta <- m$mle$theta
  }
  
  formula<-m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  aInd <- NULL
  nbAnimals <- length(unique(data$ID))
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    g0covs <- model.matrix(recharge$g0,data[aInd,])
    recovs <- model.matrix(recharge$theta,data)
    nbRecovs <- ncol(recovs)-1
    data$recharge<-rep(0,nrow(data))
    for(i in 1:nbAnimals){
      idInd <- which(data$ID==unique(data$ID)[i])
      if(nbRecovs){
        g <- g0 %*% t(g0covs[i,,drop=FALSE])
        data$recharge[idInd] <- cumsum(c(g,theta%*%t(recovs[idInd[-length(idInd)],])))
      }
    }
    newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
  }
  
  nbCovs <- ncol(model.matrix(newformula,data))-1 # substract intercept column
  
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(!isFALSE(m$conditions$circularAngleMean[[i]])) {
      meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  
  consensus <- vector('list',length(distnames))
  names(consensus) <- distnames
  for(i in distnames){
    consensus[[i]] <- (dist[[i]]=="vmConsensus")
  }
  dist <- lapply(dist,function(x) gsub("Consensus","",x))

  par <- w2n(m$mod$estimate,m$conditions$bounds,lapply(m$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,consensus,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,nbObs,dist,m$conditions$Bndind,nc,meanind,m$covsDelta,m$conditions$workBounds)
  
  Fun <- lapply(dist,function(x) paste("d",x,sep=""))
  

  probs <- matrix(1,nrow=nbObs,ncol=nbStates)
  
  for(i in distnames){
  
    genInd <- which(!is.na(data[[i]]))
    sp <- par[[i]]
  
    for(state in 1:nbStates) {
      genPar <- sp
      genProb <- rep(1,nbObs)
      genFun <- Fun[[i]]
      
      # Constitute the lists of state-dependent parameters for the step and angle
      genArgs <- list(data[[i]][genInd])
      
      zeromass <- 0
      onemass <- 0
      if(zeroInflation[[i]] | oneInflation[[i]]) {
        if(zeroInflation[[i]]) zeromass <- genPar[nrow(genPar)-nbStates*oneInflation[[i]]-nbStates+state,genInd]
        if(oneInflation[[i]]) onemass <- genPar[nrow(genPar)-nbStates+state,genInd]
        genPar <- genPar[-(nrow(genPar)-(nbStates*(zeroInflation[[i]]+oneInflation[[i]])-1):0),]
      }
  
      for(j in 1:(nrow(genPar)/nbStates))
        genArgs[[j+1]] <- genPar[(j-1)*nbStates+state,genInd]
  
      # conversion between mean/sd and shape/scale if necessary
      if(dist[[i]]=="gamma") {
        shape <- genArgs[[2]]^2/genArgs[[3]]^2
        scale <- genArgs[[3]]^2/genArgs[[2]]
        genArgs[[2]] <- shape
        genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
      }
      if(zeroInflation[[i]] | oneInflation[[i]]) {
        if(zeroInflation[[i]] & !oneInflation[[i]]){
          genProb[genInd] <- ifelse(data[[i]][genInd]==0,
                                      zeromass, # if gen==0
                                      (1-zeromass)*do.call(genFun,genArgs)) # if gen != 0
        } else if(oneInflation[[i]] & !zeroInflation[[i]]){
          genProb[genInd] <- ifelse(data[[i]][genInd]==1,
                                    onemass, # if gen==0
                                    (1-onemass)*do.call(genFun,genArgs)) # if gen != 1          
        } else {
          genProb[genInd][data[[i]][genInd]==0] <- zeromass[data[[i]][genInd]==0]
          genProb[genInd][data[[i]][genInd]==1] <- onemass[data[[i]][genInd]==1]
          genProb[genInd][data[[i]][genInd]>0 & data[[i]][genInd]<1] <- (1.-zeromass[data[[i]][genInd]>0 & data[[i]][genInd]<1]-onemass[data[[i]][genInd]>0 & data[[i]][genInd]<1]) * do.call(genFun,genArgs)[data[[i]][genInd]>0 & data[[i]][genInd]<1] # if gen !=0 and gen!=1
        }
      }
      else genProb[genInd] <- do.call(genFun,genArgs)
  
      probs[,state] <- probs[,state]*genProb;
    }
  }
  if(!is.null(m$knownStates)) {
    for (i in which(!is.na(m$knownStates))) {
      prob <- probs[i, m$knownStates[i]]
      probs[i, ] <- 0
      probs[i, m$knownStates[i]] <- prob
    }
  }
  return(probs)
}


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
#' 
#' @importFrom extraDistr dcat
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
      
      m$conditions$workBounds[[i]]<-matrix(c(-Inf,Inf),nrow(m$conditions$workBounds[[i]]),2,byrow=TRUE)
    }
    
    Par<-lapply(Par[distnames],function(x) c(t(x)))
    beta <- m$Par$beta$beta$est
    pie <- m$Par$real$pi$est
    delta <- m$Par$real$delta$est
    if(!is.null(beta)) m$conditions$workBounds$beta<-matrix(c(-Inf,Inf),length(beta),2,byrow=TRUE)
    if(!is.null(pie)) m$conditions$workBounds$pi <- matrix(c(-Inf,Inf),length(m$Par$beta$pi$est),2,byrow=TRUE)
    if(!is.null(m$Par$beta$delta$est)) m$conditions$workBounds$delta<-matrix(c(-Inf,Inf),length(m$Par$beta$delta$est),2,byrow=TRUE)
    
    g0 <- c(m$Par$beta$g0$est)
    theta <- c(m$Par$beta$theta$est)
    if(!is.null(g0)) m$conditions$workBounds$g0<-matrix(c(-Inf,Inf),length(g0),2,byrow=TRUE)
    if(!is.null(theta)) m$conditions$workBounds$theta<-matrix(c(-Inf,Inf),length(theta),2,byrow=TRUE)
    
    inputs <- checkInputs(nbStates,dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$stateNames)
    p <- inputs$p
    DMinputs<-getDM(data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
    m$conditions$fullDM <- DMinputs$fullDM
    m$mod$estimate <- n2w(Par,p$bounds,list(beta=beta,pi=m$Par$beta$pi$est,g0=g0,theta=theta),m$Par$beta$delta$est,nbStates,inputs$estAngleMean,inputs$DM,p$Bndind,inputs$dist)
  } else {
    beta <- m$mle$beta
    pie <- m$mle$pi
    delta <- m$mle$delta
    g0 <- m$mle$g0
    theta <- m$mle$theta
  }
  
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,data,par=m$mle)
  data <- cbind(data,reForm$newdata)
  nbCovs <- reForm$nbCovs
  
  ncmean <- get_ncmean(distnames,m$conditions$fullDM,m$conditions$circularAngleMean,nbStates)
  nc <- ncmean$nc
  meanind <- ncmean$meanind
  
  consensus <- vector('list',length(distnames))
  names(consensus) <- distnames
  for(i in distnames){
    consensus[[i]] <- (dist[[i]]=="vmConsensus")
  }
  dist <- lapply(dist,function(x) gsub("Consensus","",x))
  dist <- lapply(dist,function(x) ifelse(grepl("cat",x),"cat",x))

  par <- w2n(m$mod$estimate,m$conditions$bounds,lapply(m$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,consensus,m$conditions$stationary,m$conditions$fullDM,m$conditions$DMind,nbObs,dist,m$conditions$Bndind,nc,meanind,m$covsDelta,m$conditions$workBounds,m$covsPi)
  
  Fun <- lapply(dist,function(x) paste("d",x,sep=""))
  

  probs <- matrix(1,nrow=nbObs,ncol=nbStates)
  
  for(i in distnames){
    
    if(dist[[i]] %in% mvndists){
      if(dist[[i]]=="mvnorm2" || dist[[i]]=="rw_mvnorm2")
        genData <- c(data[[paste0(i,".x")]],data[[paste0(i,".y")]])
      else if(dist[[i]]=="mvnorm3" || dist[[i]]=="rw_mvnorm3")
        genData <- c(data[[paste0(i,".x")]],data[[paste0(i,".y")]],data[[paste0(i,".z")]])
    } else {
      genData <- data[[i]]
    }
  
    genInd <- which(!is.na(genData[1:nbObs]))
    sp <- par[[i]]
  
    for(state in 1:nbStates) {
      genPar <- sp
      genProb <- rep(1,nbObs)
      genFun <- Fun[[i]]
      
      # Constitute the lists of state-dependent parameters for the step and angle
      genArgs <- list(genData[which(!is.na(genData))])
      
      zeromass <- 0
      onemass <- 0
      if(zeroInflation[[i]] | oneInflation[[i]]) {
        if(zeroInflation[[i]]) zeromass <- genPar[nrow(genPar)-nbStates*oneInflation[[i]]-nbStates+state,genInd]
        if(oneInflation[[i]]) onemass <- genPar[nrow(genPar)-nbStates+state,genInd]
        genPar <- genPar[-(nrow(genPar)-(nbStates*(zeroInflation[[i]]+oneInflation[[i]])-1):0),]
      }
  
      if(dist[[i]] %in% mvndists){
        if(dist[[i]]=="mvnorm2" || dist[[i]]=="rw_mvnorm2"){
          genArgs[[2]] <- rbind(genPar[state,genInd],
                                genPar[nbStates+state,genInd])
          genArgs[[3]] <- rbind(genPar[nbStates*2+state,genInd], #x
                                genPar[nbStates*3+state,genInd], #xy
                                genPar[nbStates*3+state,genInd], #xy
                                genPar[nbStates*4+state,genInd]) #y
        } else if(dist[[i]]=="mvnorm3" || dist[[i]]=="rw_mvnorm3"){
          genArgs[[2]] <- rbind(genPar[state,genInd],
                                genPar[nbStates+state,genInd],
                                genPar[2*nbStates+state,genInd])
          genArgs[[3]] <- rbind(genPar[nbStates*3+state,genInd], #x
                                genPar[nbStates*4+state,genInd], #xy
                                genPar[nbStates*5+state,genInd], #xz
                                genPar[nbStates*4+state,genInd], #xy
                                genPar[nbStates*6+state,genInd], #y
                                genPar[nbStates*7+state,genInd], #yz
                                genPar[nbStates*5+state,genInd], #xz
                                genPar[nbStates*7+state,genInd], #yz
                                genPar[nbStates*8+state,genInd]) #z          
        }
      } else if(dist[[i]]=="cat"){
        dimCat <- as.numeric(gsub("cat","",m$conditions$dist[[i]]))
        genArgs[[2]] <- t(genPar[seq(state,dimCat*nbStates,nbStates),genInd])
      } else {
        for(j in 1:(nrow(genPar)/nbStates))
          genArgs[[j+1]] <- genPar[(j-1)*nbStates+state,genInd]
      }
      
      # conversion between mean/sd and shape/scale if necessary
      if(dist[[i]]=="gamma") {
        shape <- genArgs[[2]]^2/genArgs[[3]]^2
        scale <- genArgs[[3]]^2/genArgs[[2]]
        genArgs[[2]] <- shape
        genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
      }
      if(zeroInflation[[i]] | oneInflation[[i]]) {
        if(zeroInflation[[i]] & !oneInflation[[i]]){
          genProb[genInd] <- ifelse(genData[genInd]==0,
                                      zeromass, # if gen==0
                                      (1-zeromass)*do.call(genFun,genArgs)) # if gen != 0
        } else if(oneInflation[[i]] & !zeroInflation[[i]]){
          genProb[genInd] <- ifelse(genData[genInd]==1,
                                    onemass, # if gen==0
                                    (1-onemass)*do.call(genFun,genArgs)) # if gen != 1          
        } else {
          genProb[genInd][genData[genInd]==0] <- zeromass[genData[genInd]==0]
          genProb[genInd][genData[genInd]==1] <- onemass[genData[genInd]==1]
          genProb[genInd][genData[genInd]>0 & genData[genInd]<1] <- (1.-zeromass[genData[genInd]>0 & genData[genInd]<1]-onemass[genData[genInd]>0 & genData[genInd]<1]) * do.call(genFun,genArgs)[genData[genInd]>0 & genData[genInd]<1] # if gen !=0 and gen!=1
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

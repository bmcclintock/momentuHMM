
#' Pseudo-residuals
#'
#' The pseudo-residuals of momentuHMM models, as described in Zucchini and McDonad (2009).
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, \code{\link{HMMfits}}, or \code{\link{miSum}} object.
#' @param ncores number of cores to use for parallel processing
#'
#' @return If \code{m} is a \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object, a list of pseudo-residuals for each data stream (e.g., 'stepRes', 'angleRes') is returned. 
#' If \code{m} is a list of \code{\link{momentuHMM}} objects, then a list of length \code{length(m)} is returned where each element is a list of pseudo-residuals for each data stream.
#'
#' @details If some turning angles in the data are equal to pi, the corresponding pseudo-residuals
#' will not be included. Indeed, given that the turning angles are defined on (-pi,pi], an angle of pi
#' results in a pseudo-residual of +Inf (check Section 6.2 of reference for more information on the
#' computation of pseudo-residuals).
#' 
#' A continuity adjustment (adapted from Harte 2017) is made for discrete probability distributions. When
#' the data are near the boundary (e.g. 0 for ``pois''; 0 and 1 for ``bern''), then the pseudo residuals can
#' be a poor indicator of lack of fit.
#' 
#' For multiple imputation analyses, if \code{m} is a \code{\link{miHMM}} object or a list of \code{\link{momentuHMM}} objects, then
#' the pseudo-residuals are individually calculated for each model fit. Note that pseudo-residuals for \code{\link{miSum}} objects (as returned by \code{\link{MIpool}}) are based on pooled parameter 
#' estimates and the means of the data values across all imputations (and therefore may not be particularly meaningful).
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' res <- pseudoRes(m)
#' qqnorm(res$stepRes)
#' qqnorm(res$angleRes)
#'
#' @references
#' Harte, D. 2017. HiddenMarkov: Hidden Markov Models. R package version 1.8-8.
#'
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export
#' @importFrom stats integrate qnorm qchisq mahalanobis
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom extraDistr pcat

pseudoRes <- function(m, ncores = 1)
{

  m <- delta_bc(m)
  
  if(!is.momentuHMM(m) & !is.miSum(m)){
    if(!is.miHMM(m) & !is.HMMfits(m)) stop("'m' must be a momentuHMM, HMMfits, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
    else {
      if(is.miHMM(m)) m <- m$HMMfits
      registerDoParallel(cores=ncores)
      withCallingHandlers(genRes <- foreach(i=which(unlist(lapply(m,is.momentuHMM)))) %dorng% {
        pseudoRes(m[[i]])
      },warning=muffleRNGwarning)
      stopImplicitCluster()
      return(genRes)
    }
  }

  data <- m$data
  nbObs <- nrow(data)
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))
  dist <- m$conditions$dist
  distnames <- names(dist)
  
  if(is.miSum(m)){
    warning('pseudo-residuals are based on pooled parameter estimates and mean data values across multiple imputations...')
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
    
    m$mle <- list(g0=g0,theta=theta)
    
    inputs <- checkInputs(nbStates,dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
    p <- inputs$p
    DMinputs<-getDM(data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
    m$conditions$fullDM <- DMinputs$fullDM
    m$mod$estimate <- n2w(Par,p$bounds,list(beta=beta,pi=m$Par$beta$pi$est,g0=g0,theta=theta),m$Par$beta$delta$est,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind,inputs$dist)
  } else {
    beta <- m$mle$beta
    pie <- m$mle$pi
    delta <- m$mle$delta
    g0 <- m$mle$g0
    theta <- m$mle$theta
  }
  
  consensus <- vector('list',length(distnames))
  names(consensus) <- distnames
  
  for(i in distnames){
    if(dist[[i]]=="vmConsensus"){
      consensus[[i]] <- TRUE
      dist[[i]] <- gsub("Consensus","",dist[[i]])
    } else consensus[[i]] <- FALSE
  }
  dist <- lapply(dist,function(x) ifelse(grepl("cat",x),"cat",x))
  
  Fun <- lapply(dist,function(x) paste("p",x,sep=""))
  for(j in which(dist %in% angledists)){
    Fun[[j]] <- paste0("d",dist[[j]])
    if(length(which(data[[distnames[j]]]==pi))>0)
      message("Note: Some ",distnames[j],"s are equal to pi, and the corresponding pseudo-residuals are not included")
  }

  # forward log-probabilities
  la <- logAlpha(m)
  
  # identify covariates
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,data,par=list(g0=g0,theta=theta))
  recharge <- reForm$recharge
  hierRecharge <- reForm$hierRecharge
  newformula <- reForm$newformula
  if(!is.null(recharge)) data[colnames(reForm$newdata)] <- reForm$newdata
  covs <- reForm$covs
  nbCovs <- reForm$nbCovs
  aInd <- reForm$aInd
  
  mixtures <- m$conditions$mixtures
  if(mixtures==1) pie <- matrix(1,nbAnimals,1)
  
  ncmean <- get_ncmean(distnames,m$conditions$fullDM,m$conditions$circularAngleMean,nbStates)
  nc <- ncmean$nc
  meanind <- ncmean$meanind
  
  par <- w2n(m$mod$estimate,m$conditions$bounds,lapply(m$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,consensus,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,nbObs,dist,m$conditions$Bndind,nc,meanind,m$covsDelta,m$conditions$workBounds,m$covsPi)
  
  trMat <- list()
  for(mix in 1:mixtures){
    if(nbStates>1){
      if(is.null(recharge)){
        trMat[[mix]] <- trMatrix_rcpp(nbStates,beta[(mix-1)*(nbCovs+1)+1:(nbCovs+1),,drop=FALSE],as.matrix(covs),m$conditions$betaRef)
      } else {
        gamInd<-(length(m$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)*mixtures+1):(length(m$mod$estimate))-(ncol(m$covsPi)*(mixtures-1))-ifelse(reForm$nbRecovs,reForm$nbRecovs+1+reForm$nbG0covs+1,0)-ncol(m$covsDelta)*(nbStates-1)*(!m$conditions$stationary)*mixtures
        trMat[[mix]] <- array(unlist(lapply(split(data,1:nrow(data)),function(x) tryCatch(get_gamma_recharge(m$mod$estimate[c(gamInd[unique(c(m$conditions$betaCons))],length(m$mod$estimate)-reForm$nbRecovs:0)],covs=x,formula=newformula,hierRecharge=hierRecharge,nbStates=nbStates,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=rbind(m$conditions$workBounds$beta,m$conditions$workBounds$theta),mixture=mix),error=function(e) NA))),dim=c(nbStates,nbStates,nrow(data))) 
      }
    } else trMat[[mix]] <- array(1,dim=c(1,1,nbObs))
  }
  
  genRes <- list()
  for(j in distnames){
    genRes[[paste0(j,"Res")]] <- rep(0,nbObs)
    pgenMat <- pgenMat2 <- matrix(NA,nbObs,nbStates)
    sp <- par[[j]]
    
    if(dist[[j]] %in% mvndists){
      if(dist[[j]]=="mvnorm2" || dist[[j]]=="rw_mvnorm2"){
        genData <- c(data[[paste0(j,".x")]],data[[paste0(j,".y")]])
        if(dist[[j]]=="mvnorm2") ndim <- as.numeric(gsub("mvnorm","",dist[[j]]))
        else ndim <- as.numeric(gsub("rw_mvnorm","",dist[[j]]))
      } else if(dist[[j]]=="mvnorm3" || dist[[j]]=="rw_mvnorm3"){
        genData <- c(data[[paste0(j,".x")]],data[[paste0(j,".y")]],data[[paste0(j,".z")]])
        if(dist[[j]]=="mvnorm3") ndim <- as.numeric(gsub("mvnorm","",dist[[j]]))
        else ndim <- as.numeric(gsub("rw_mvnorm","",dist[[j]]))
      }
      
      # define function for calculating chi square probs based on mahalanobis distance
      d2<-function(q,mean,sigma){
        stats::pchisq(stats::mahalanobis(q,mean,matrix(sigma,length(mean),length(mean))),df=ndim)
      }
      
    } else {
      genData <- data[[j]]
    }
    genInd <- which(!is.na(genData[1:nbObs]))
    zeroInflation <- m$conditions$zeroInflation[[j]]
    oneInflation <- m$conditions$oneInflation[[j]]
  
    for(state in 1:nbStates) {
      
      genPar <- sp
      
      if(!(dist[[j]] %in% angledists)){
        
        genArgs <- list(genData[which(!is.na(genData))])
        
        zeromass <- 0
        onemass <- 0
        if(zeroInflation | oneInflation) {
          if(zeroInflation) zeromass <- genPar[nrow(genPar)-nbStates*oneInflation-nbStates+state,genInd]
          if(oneInflation) onemass <- genPar[nrow(genPar)-nbStates+state,genInd]
          genPar <- genPar[-(nrow(genPar)-(nbStates*(zeroInflation+oneInflation)-1):0),]
        }
        
        if(dist[[j]] %in% mvndists){
          if(dist[[j]]=="mvnorm2" || dist[[j]]=="rw_mvnorm2"){
            genArgs[[1]] <- as.list(as.data.frame(matrix(genArgs[[1]],nrow=ndim,byrow=TRUE)))
            genArgs[[2]] <- as.list(as.data.frame(rbind(genPar[state,genInd],
                                  genPar[nbStates+state,genInd])))
            genArgs[[3]] <- as.list(as.data.frame(
                                  rbind(genPar[nbStates*2+state,genInd], #x
                                        genPar[nbStates*3+state,genInd], #xy
                                        genPar[nbStates*3+state,genInd], #xy
                                        genPar[nbStates*4+state,genInd]))) #y
          } else if(dist[[j]]=="mvnorm3" || dist[[j]]=="rw_mvnorm3"){
            genArgs[[1]] <- as.list(as.data.frame(matrix(genArgs[[1]],nrow=ndim,byrow=TRUE)))
            genArgs[[2]] <- as.list(as.data.frame(rbind(genPar[state,genInd],
                                  genPar[nbStates+state,genInd],
                                  genPar[2*nbStates+state,genInd])))
            genArgs[[3]] <- as.list(as.data.frame(
                                  rbind(genPar[nbStates*3+state,genInd], #x
                                        genPar[nbStates*4+state,genInd], #xy
                                        genPar[nbStates*5+state,genInd], #xz
                                        genPar[nbStates*4+state,genInd], #xy
                                        genPar[nbStates*6+state,genInd], #y
                                        genPar[nbStates*7+state,genInd], #yz
                                        genPar[nbStates*5+state,genInd], #xz
                                        genPar[nbStates*7+state,genInd], #yz
                                        genPar[nbStates*8+state,genInd]))) #z          
          }
        } else if(dist[[j]]=="cat"){
          dimCat <- as.numeric(gsub("cat","",m$conditions$dist[[j]]))
          genArgs[[2]] <- t(genPar[seq(state,dimCat*nbStates,nbStates),genInd])
        } else {
          for(k in 1:(nrow(genPar)/nbStates))
            genArgs[[k+1]] <- genPar[(k-1)*nbStates+state,genInd]
        }
        
        if(dist[[j]]=="gamma") {
          shape <- genArgs[[2]]^2/genArgs[[3]]^2
          scale <- genArgs[[3]]^2/genArgs[[2]]
          genArgs[[2]] <- shape
          genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        
        if(zeroInflation | oneInflation) {
          if(zeroInflation & !oneInflation){
            pgenMat[genInd,state] <- ifelse(genData[genInd]==0,
                                      zeromass, # if gen==0
                                      zeromass + (1-zeromass)*do.call(Fun[[j]],genArgs)) # if gen != 0
          } else if(oneInflation & !zeroInflation){
            pgenMat[genInd,state] <- ifelse(genData[genInd]==1,
                                      onemass, # if gen==1
                                      onemass + (1-onemass)*do.call(Fun[[j]],genArgs)) # if gen != 1           
          } else {
            pgenMat[genInd,state][genData[genInd]==0] <- zeromass[genData[genInd]==0] # if gen==0
            pgenMat[genInd,state][genData[genInd]==1] <- onemass[genData[genInd]==1]  # if gen==1
            pgenMat[genInd,state][genData[genInd]>0 & genData[genInd]<1] <- zeromass[genData[genInd]>0 & genData[genInd]<1] + onemass[genData[genInd]>0 & genData[genInd]<1] + (1.-zeromass[genData[genInd]>0 & genData[genInd]<1]-onemass[genData[genInd]>0 & genData[genInd]<1]) * do.call(Fun[[j]],genArgs)[genData[genInd]>0 & genData[genInd]<1] # if gen !=0 and gen!=1
          }
        }
        else {
          if(dist[[j]] %in% mvndists){
            names(genArgs) <- c("q","mean","sigma")
            pgenMat[genInd,state] <- mapply(d2,q=genArgs$q,mean=genArgs$mean,sigma=genArgs$sigma)
          } else pgenMat[genInd,state] <- do.call(Fun[[j]],genArgs)
          if(dist[[j]] %in% integerdists){
            genArgs[[1]] <- genArgs[[1]] - 1
            pgenMat2[genInd,state] <- do.call(Fun[[j]],genArgs)
          }
        }
        
      } else {
        
        genpiInd <- which(genData!=pi & !is.na(genData))
        
        genArgs <- list(Fun[[j]],-pi,genData[1]) # to pass to function "integrate" below
  
        for(i in genpiInd){
          genArgs[[3]]<-genData[i]
          for(k in 1:(nrow(genPar)/nbStates))
            genArgs[[k+3]] <- genPar[(k-1)*nbStates+state,i]
          
          pgenMat[i,state] <- do.call(integrate,genArgs)$value
        }
      }
    }
    
    k <- 1
    for(i in 1:nbObs) {
      if(any(i==aInd)) {
        iPi <- pie[k,]
        kInd <- k
        k <- k + 1
      }
      if(!is.na(genData[i])){
        for(mix in 1:mixtures){
          if(any(i==aInd)){
            #iPi <- pie[k,]
            if(dist[[j]] %in% integerdists)
              genRes[[paste0(j,"Res")]][i] <- genRes[[paste0(j,"Res")]][i] + ((delta[(mix-1)*nbAnimals+kInd,]%*%trMat[[mix]][,,i])%*%pgenMat[i,] + (delta[(mix-1)*nbAnimals+kInd,]%*%trMat[[mix]][,,i])%*%pgenMat2[i,])/2 * iPi[mix]
            else
              genRes[[paste0(j,"Res")]][i] <- genRes[[paste0(j,"Res")]][i] + (delta[(mix-1)*nbAnimals+kInd,]%*%trMat[[mix]][,,i])%*%pgenMat[i,] * iPi[mix]
          } else {
            gamma <- trMat[[mix]][,,i]
            #c <- max(la[i-1,]) # cancels below ; prevents numerical errors
            #a <- exp(la[i-1,]-c)
            c <- max(la[[mix]][i-1,])
            a <- exp(la[[mix]][i-1,]-c)
            
            if(dist[[j]] %in% integerdists)
              genRes[[paste0(j,"Res")]][i] <- genRes[[paste0(j,"Res")]][i] + (t(a)%*%(gamma/sum(a))%*%pgenMat[i,] + t(a)%*%(gamma/sum(a))%*%pgenMat2[i,])/2 * iPi[mix]
            else
              genRes[[paste0(j,"Res")]][i] <- genRes[[paste0(j,"Res")]][i] + t(a)%*%(gamma/sum(a))%*%pgenMat[i,] * iPi[mix]
          }
        }
        if(dist[[j]] %in% mvndists){
          genRes[[paste0(j,"Res")]][i] <- stats::qchisq(genRes[[paste0(j,"Res")]][i],df=ndim)
        } else genRes[[paste0(j,"Res")]][i] <- stats::qnorm(genRes[[paste0(j,"Res")]][i])
      } else genRes[[paste0(j,"Res")]][i] <- NA
    }
  }

  return(genRes)
}

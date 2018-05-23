
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
#' @importFrom stats integrate qnorm
#' @importFrom LaplacesDemon pbern
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%

pseudoRes <- function(m, ncores = 1)
{

  m <- delta_bc(m)
  
  if(!is.momentuHMM(m) & !is.miSum(m)){
    if(!is.miHMM(m) & !is.HMMfits(m)) stop("'m' must be a momentuHMM, HMMfits, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
    else {
      if(is.miHMM(m)) m <- m$HMMfits
      registerDoParallel(cores=ncores)
      genRes <- foreach(i=which(unlist(lapply(m,is.momentuHMM)))) %dopar% {
        pseudoRes(m[[i]])
      }
      stopImplicitCluster()
      return(genRes)
    }
  }

  data <- m$data
  nbObs <- nrow(data)
  nbStates <- length(m$stateNames)
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
    if(!is.null(m$mle$beta)) m$conditions$workBounds$beta<-matrix(c(-Inf,Inf),length(m$mle$beta),2,byrow=TRUE)
    if(!is.null(m$Par$beta$delta$est)) m$conditions$workBounds$delta<-matrix(c(-Inf,Inf),length(m$Par$beta$delta$est),2,byrow=TRUE)
    
    Par<-lapply(Par,function(x) c(t(x)))
    Par<-Par[distnames]
    beta <- m$Par$beta$beta$est
    delta <- m$Par$real$delta$est
    inputs <- checkInputs(nbStates,dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
    p <- inputs$p
    DMinputs<-getDM(data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
    m$conditions$fullDM <- DMinputs$fullDM
    m$mod$estimate <- n2w(Par,p$bounds,beta,m$Par$beta$delta$est,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
  }
  
  consensus <- vector('list',length(distnames))
  names(consensus) <- distnames
  
  for(i in distnames){
    if(dist[[i]]=="vmConsensus"){
      consensus[[i]] <- TRUE
      dist[[i]] <- gsub("Consensus","",dist[[i]])
    } else consensus[[i]] <- FALSE
  }
  
  Fun <- lapply(dist,function(x) paste("p",x,sep=""))
  for(j in which(dist %in% angledists)){
    Fun[[j]] <- paste0("d",dist[[j]])
    if(length(which(data[[distnames[j]]]==pi))>0)
      message("Note: Some ",distnames[j],"s are equal to pi, and the corresponding pseudo-residuals are not included")
  }

  # forward log-probabilities
  la <- logAlpha(m)
  
  # identify covariates
  formula<-m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  
  covs <- model.matrix(newformula,data)
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  aInd <- NULL
  for(i in 1:length(unique(m$data$ID)))
    aInd <- c(aInd,which(m$data$ID==unique(m$data$ID)[i])[1])
  
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) {
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
  
  par <- w2n(m$mod$estimate,m$conditions$bounds,lapply(m$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,consensus,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,nbObs,dist,m$conditions$Bndind,nc,meanind,m$covsDelta,m$conditions$workBounds)
  
  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  genRes <- list()
  for(j in distnames){
    genRes[[paste0(j,"Res")]] <- rep(NA,nbObs)
    pgenMat <- pgenMat2 <- matrix(NA,nbObs,nbStates)
    sp <- par[[j]]
    genInd <- which(!is.na(data[[j]]))
    zeroInflation <- m$conditions$zeroInflation[[j]]
    oneInflation <- m$conditions$oneInflation[[j]]
  
    for(state in 1:nbStates) {
      
      genPar <- sp
      
      if(!(dist[[j]] %in% angledists)){
        
        genArgs <- list(data[[j]][genInd])
        
        zeromass <- 0
        onemass <- 0
        if(zeroInflation | oneInflation) {
          if(zeroInflation) zeromass <- genPar[nrow(genPar)-nbStates*oneInflation-nbStates+state,genInd]
          if(oneInflation) onemass <- genPar[nrow(genPar)-nbStates+state,genInd]
          genPar <- genPar[-(nrow(genPar)-(nbStates*(zeroInflation+oneInflation)-1):0),]
        }
        for(k in 1:(nrow(genPar)/nbStates))
          genArgs[[k+1]] <- genPar[(k-1)*nbStates+state,genInd]
        
        if(dist[[j]]=="gamma") {
          shape <- genArgs[[2]]^2/genArgs[[3]]^2
          scale <- genArgs[[3]]^2/genArgs[[2]]
          genArgs[[2]] <- shape
          genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        
        if(zeroInflation | oneInflation) {
          if(zeroInflation & !oneInflation){
            pgenMat[genInd,state] <- ifelse(data[[j]][genInd]==0,
                                      zeromass, # if gen==0
                                      zeromass + (1-zeromass)*do.call(Fun[[j]],genArgs)) # if gen != 0
          } else if(oneInflation & !zeroInflation){
            pgenMat[genInd,state] <- ifelse(data[[j]][genInd]==1,
                                      onemass, # if gen==1
                                      onemass + (1-onemass)*do.call(Fun[[j]],genArgs)) # if gen != 1           
          } else {
            pgenMat[genInd,state][data[[j]][genInd]==0] <- zeromass[data[[j]][genInd]==0] # if gen==0
            pgenMat[genInd,state][data[[j]][genInd]==1] <- onemass[data[[j]][genInd]==1]  # if gen==1
            pgenMat[genInd,state][data[[j]][genInd]>0 & data[[j]][genInd]<1] <- zeromass[data[[j]][genInd]>0 & data[[j]][genInd]<1] + onemass[data[[j]][genInd]>0 & data[[j]][genInd]<1] + (1.-zeromass[data[[j]][genInd]>0 & data[[j]][genInd]<1]-onemass[data[[j]][genInd]>0 & data[[j]][genInd]<1]) * do.call(Fun[[j]],genArgs)[data[[j]][genInd]>0 & data[[j]][genInd]<1] # if gen !=0 and gen!=1
          }
        }
        else {
          pgenMat[genInd,state] <- do.call(Fun[[j]],genArgs)
          if(dist[[j]] %in% integerdists){
            genArgs[[1]] <- genArgs[[1]] - 1
            pgenMat2[genInd,state] <- do.call(Fun[[j]],genArgs)
          }
        }
        
        #pgenMat[genInd,state] <- zeromass+(1-zeromass)*do.call(Fun[[j]],genArgs)
        #for(i in 1:nbObs) {
        #  if(!is.na(data[[j]][i])) {
        #    genArgs[[1]] <- data[[j]][i]
        #    pgenMat[i,state] <- zeromass+(1-zeromass)*do.call(Fun[[j]],genArgs)
        #  }
        #}
      } else {
        
        genpiInd <- which(data[[j]]!=pi & !is.na(data[[j]]))
        
        genArgs <- list(Fun[[j]],-pi,data[[j]][1]) # to pass to function "integrate" below
  
        for(i in genpiInd){
          genArgs[[3]]<-data[[j]][i]
          for(k in 1:(nrow(genPar)/nbStates))
            genArgs[[k+3]] <- genPar[(k-1)*nbStates+state,i]
          
          pgenMat[i,state] <- do.call(integrate,genArgs)$value
        }
        #for(i in 1:nbObs) {
        #  if(!is.na(data[[j]][i])) {
        #    # angle==pi => residual=Inf
        #    if(data[[j]][i]!=pi) {
        #      genArgs[[3]] <- data[[j]][i]
        #      pgenMat[i,state] <- do.call(integrate,genArgs)$value
        #    }
        #  }
        #}
      }
    }
  
    k <- 1
    for(i in 1:nbObs) {
      if(!is.na(data[[j]][i])){
        if(any(i==aInd)){
          genRes[[paste0(j,"Res")]][i] <- (delta[k,]%*%trMat[,,i])%*%pgenMat[i,]
          if(dist[[j]] %in% integerdists)
            genRes[[paste0(j,"Res")]][i] <- (genRes[[paste0(j,"Res")]][i] + (delta[k,]%*%trMat[,,i])%*%pgenMat2[i,])/2
          genRes[[paste0(j,"Res")]][i] <- qnorm(genRes[[paste0(j,"Res")]][i])
          k <- k + 1
        } else {
          gamma <- trMat[,,i]
          c <- max(la[i-1,]) # cancels below ; prevents numerical errors
          a <- exp(la[i-1,]-c)
          
          genRes[[paste0(j,"Res")]][i] <- t(a)%*%(gamma/sum(a))%*%pgenMat[i,]
          if(dist[[j]] %in% integerdists)
            genRes[[paste0(j,"Res")]][i] <- (genRes[[paste0(j,"Res")]][i] + t(a)%*%(gamma/sum(a))%*%pgenMat2[i,])/2
          genRes[[paste0(j,"Res")]][i] <- qnorm(genRes[[paste0(j,"Res")]][i])
        }
      }
    }
  }

  return(genRes)
}

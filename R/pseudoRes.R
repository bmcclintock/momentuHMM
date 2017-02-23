
#' Pseudo-residuals
#'
#' The pseudo-residuals of a momentuHMM model, as described in Zucchini and McDonad (2009).
#'
#' @param m A \code{momentuHMM} object.
#'
#' @return A list of:
#' \item{stepRes}{The pseudo-residuals for the step lengths}
#' \item{angleRes}{The pseudo-residuals for the turning angles}
#'
#' @details If some turning angles in the data are equal to pi, the corresponding pseudo-residuals
#' will not be included. Indeed, given that the turning angles are defined on (-pi,pi], an angle of pi
#' results in a pseudo-residual of +Inf (check Section 6.2 of reference for more information on the
#' computation of pseudo-residuals).
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' res <- pseudoRes(m)
#' qqnorm(res$stepRes)
#' qqnorm(res$angleRes)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export
#' @importFrom stats integrate qnorm

pseudoRes <- function(m)
{
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(m)) m <- m$miSum

  data <- m$data
  nbObs <- nrow(data)
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  
  if(is.miSum(m)){
    warning('pseudo-residuals are based on pooled parameter estimates and mean covariate values across multiple imputations...')
    Par <- lapply(m$Par$real,function(x) x$est)
    for(i in distnames){
      if(!is.null(m$conditions$DM[[i]]))
        Par[[i]] <- m$Par$beta[[i]]$est
      else if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]])
        Par[[i]] <- Par[[i]][-1,]
    }
    Par<-lapply(Par,function(x) c(t(x)))
    Par<-Par[distnames]
    beta <- m$Par$beta$beta$est
    delta <- m$Par$real$delta$est
    inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
    p <- inputs$p
    DMinputs<-getDM(data,inputs$DM,m$conditions$dist,nbStates,p$parNames,p$bounds,Par,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$circularAngleMean[i])
    m$conditions$fullDM <- DMinputs$fullDM
    m$mod$estimate <- n2w(Par,p$bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
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
  covs <- model.matrix(m$conditions$formula,data)
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  par <- w2n(m$mod$estimate,m$conditions$bounds,lapply(m$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,nbObs,dist,m$conditions$Bndind)
  
  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  genRes <- list()
  for(j in distnames){
    genRes[[paste0(j,"Res")]] <- rep(NA,nbObs)
    pgenMat <- matrix(NA,nbObs,nbStates)
    sp <- par[[j]]
    genInd <- which(!is.na(data[[j]]))
  
    for(state in 1:nbStates) {
      
      genPar <- sp
      
      if(!(dist[[j]] %in% angledists)){
        
        genArgs <- list(data[[j]][genInd])
        
        if(!m$conditions$zeroInflation[[j]]) {
          zeromass <- 0
        }
        else {
          zeromass <- genPar[nrow(genPar)-nbStates+state,genInd]
          genPar <- genPar[-(nrow(genPar)-(nbStates-1):0),]
        }
        for(k in 1:(nrow(genPar)/nbStates))
          genArgs[[k+1]] <- genPar[(k-1)*nbStates+state,genInd]
        
        if(dist[[j]]=="gamma") {
          shape <- genArgs[[2]]^2/genArgs[[3]]^2
          scale <- genArgs[[3]]^2/genArgs[[2]]
          genArgs[[2]] <- shape
          genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        
        pgenMat[genInd,state] <- zeromass+(1-zeromass)*do.call(Fun[[j]],genArgs)
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
  
    if(!is.na(data[[j]][1]))
      genRes[[paste0(j,"Res")]][1] <- qnorm((delta%*%trMat[,,1])%*%pgenMat[1,])

    for(i in 2:nbObs) {
      gamma <- trMat[,,i]
      c <- max(la[i-1,]) # cancels below ; prevents numerical errors
      a <- exp(la[i-1,]-c)
  
      if(!is.na(data[[j]][i]))
        genRes[[paste0(j,"Res")]][i] <-qnorm(t(a)%*%(gamma/sum(a))%*%pgenMat[i,])
    }
  }

  return(genRes)
}

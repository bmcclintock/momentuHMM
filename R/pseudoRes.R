
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
  if(!is.momentuHMM(m) & !is.momentuHMMMI(m))
    stop("'m' must be a momentuHMM or momentuHMMMI object (as output by fitHMM or MI_summary)")

  data <- m$data
  nbObs <- nrow(data)
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  
  if(is.momentuHMMMI(m)){
    m$mle<-lapply(m$Par[distnames],function(x) x$est)
    m$mle$beta<-m$Par$beta$est
    m$mle$delta<-m$Par$delta$est
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

  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  genRes <- list()
  for(j in distnames){
    genRes[[paste0(j,"Res")]] <- rep(NA,nbObs)
    pgenMat <- matrix(NA,nbObs,nbStates)
  
    for(state in 1:nbStates) {
      # define lists of parameters
      genArgs <- list(data[[j]])
      if(!m$conditions$zeroInflation[[j]]) {
          for(k in 1:nrow(m$mle[[j]]))
              genArgs[[k+1]] <- m$mle[[j]][k,state]
  
          zeromass <- 0
      }
      else {
          for(k in 1:(nrow(m$mle[[j]])-1))
              genArgs[[k+1]] <- m$mle[[j]][k,state]
  
          zeromass <- m$mle[[j]][nrow(m$mle[[j]]),state]
      }
      if(!(dist[[j]] %in% angledists)){
        if(dist[[j]]=="gamma") {
          shape <- genArgs[[2]]^2/genArgs[[3]]^2
          scale <- genArgs[[3]]^2/genArgs[[2]]
          genArgs[[2]] <- shape
          genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        
        for(i in 1:nbObs) {
          if(!is.na(data[[j]][i])) {
            genArgs[[1]] <- data[[j]][i]
            pgenMat[i,state] <- zeromass+(1-zeromass)*do.call(Fun[[j]],genArgs)
          }
        }
      } else {
        genArgs <- list(Fun[[j]],-pi,data[[j]][1]) # to pass to function "integrate" below
        for(k in 1:nrow(m$mle[[j]]))
          genArgs[[k+3]] <- m$mle[[j]][k,state]
  
        for(i in 1:nbObs) {
          if(!is.na(data[[j]][i])) {
            # angle==pi => residual=Inf
            if(data[[j]][i]!=pi) {
              genArgs[[3]] <- data[[j]][i]
              pgenMat[i,state] <- do.call(integrate,genArgs)$value
            }
          }
        }
      }
    }
  
    if(!is.na(data[[j]][1]))
      genRes[[paste0(j,"Res")]][1] <- qnorm(t(m$mle$delta)%*%pgenMat[1,])

    for(i in 2:nbObs) {
      gamma <- trMat[,,i-1]
      c <- max(la[i-1,]) # cancels below ; prevents numerical errors
      a <- exp(la[i-1,]-c)
  
      if(!is.na(data[[j]][i]))
        genRes[[paste0(j,"Res")]][i] <-qnorm(t(a)%*%(gamma/sum(a))%*%pgenMat[i,])
    }
  }

  return(genRes)
}

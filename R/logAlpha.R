
#' Forward log-probabilities
#'
#' Used in \code{\link{stateProbs}} and \code{\link{pseudoRes}}.
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object.
#'
#' @return The matrix of forward log-probabilities.
#'
#' @examples
#' \dontrun{
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' la <- momentuHMM:::logAlpha(m)
#' }

logAlpha <- function(m)
{
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(m)) m <- m$miSum
  
  m <- delta_bc(m)
  
  if(is.miSum(m)){
    beta<-m$Par$beta$beta$est
    pie<-c(m$Par$real$pi$est)
    delta<-m$Par$real$delta$est
    g0<-m$Par$beta$g0$est
    theta<-m$Par$beta$theta$est
  } else {
    beta <- m$mle$beta
    pie <- m$mle$pi
    delta <- m$mle$delta
    g0 <- m$mle$g0
    theta <- m$mle$theta
  }
  
  nbStates <- length(m$stateNames)
  nbObs <- nrow(m$data)
  lalpha <- matrix(NA,nbObs,nbStates)

  # identify covariates
  formula<-m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  aInd <- NULL
  nbAnimals <- length(unique(m$data$ID))
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(m$data$ID==unique(m$data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    g0covs <- model.matrix(recharge$g0,m$data[aInd,])
    nbG0covs <- ncol(g0covs)-1
    recovs <- model.matrix(recharge$theta,m$data)
    nbRecovs <- ncol(recovs)-1
    m$data$recharge<-rep(0,nrow(m$data))
    for(i in 1:nbAnimals){
      idInd <- which(m$data$ID==unique(m$data$ID)[i])
      if(nbRecovs){
        g <- g0 %*% t(g0covs[i,,drop=FALSE])
        m$data$recharge[idInd] <- cumsum(c(g,theta%*%t(recovs[idInd[-length(idInd)],])))
      }
    }
    for(j in 1:nbStates){
      formulaStates[[j]] <- as.formula(paste0(Reduce( paste, deparse(formulaStates[[j]]) ),"+recharge"))
    }
    formterms <- c(formterms,"recharge")
    newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
  } else {
    nbG0covs <- 0
    nbRecovs <- 0
    g0covs <- NULL
    recovs <- NULL
  }
  
  covs <- model.matrix(newformula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  probs <- allProbs(m)
  
  mixtures <- m$conditions$mixtures
  if(mixtures==1) pie <- 1
  
  trMat <- list()
  
  foo <- matrix(0,mixtures,nbStates)
  mxlscale <- numeric(mixtures)
  
  for(mix in 1:mixtures){
    if(nbStates>1)
      trMat[[mix]] <- trMatrix_rcpp(nbStates,beta[(mix-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],as.matrix(covs),m$conditions$betaRef)
    else
      trMat[[mix]] <- array(1,dim=c(1,1,nbObs))
  }
  k <- 1
  for(i in 1:nbObs) {
    for(mix in 1:mixtures){
      gamma <- trMat[[mix]][,,i]
      if(any(i==aInd)){
        foo[mix,] <- (delta[(mix-1)*nbAnimals+k,] %*% gamma)*probs[i,]
        mxlscale[mix] <- log(pie[mix])
        if(mix==mixtures) k <- k + 1
      } else {
        foo[mix,] <- (foo[mix,] %*% gamma)*probs[i,]
      }
      mxlscale[mix] <- mxlscale[mix]+log(sum(foo[mix,]))
      foo[mix,] <- foo[mix,]/sum(foo[mix,])
    }
    A <- max(mxlscale)
    lscale <- A + log(sum(exp(mxlscale-A)))
    lalpha[i,] <- log(colSums(foo * pie))+lscale
  }
  return(lalpha)
}

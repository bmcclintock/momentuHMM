
#' Viterbi algorithm
#'
#' For a given model, reconstructs the most probable states sequence,
#' using the Viterbi algorithm.
#'
#' @param m An object \code{momentuHMM}
#'
#' @return The sequence of most probable states.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' # reconstruction of states sequence
#' states <- viterbi(m)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export

viterbi <- function(m)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  m <- delta_bc(m)
  
  data <- m$data
  nbStates <- length(m$stateNames)
  beta <- m$mle$beta
  delta <- m$mle$delta
  
  if(nbStates==1)
    stop("No states to decode (nbStates=1)")

  # identify covariates
  formula<-m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  aInd <- NULL
  nbAnimals <- length(unique(data$ID))
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    g0covs <- model.matrix(recharge$g0,data[aInd,])
    nbG0covs <- ncol(g0covs)-1
    recovs <- model.matrix(recharge$theta,data)
    nbRecovs <- ncol(recovs)-1
    data$recharge<-rep(0,nrow(data))
    for(i in 1:nbAnimals){
      idInd <- which(data$ID==unique(data$ID)[i])
      if(nbRecovs){
        g0 <- m$mle$g0 %*% t(g0covs[i,,drop=FALSE])
        theta <- m$mle$theta
        data$recharge[idInd] <- cumsum(c(g0,theta%*%t(recovs[idInd[-length(idInd)],])))
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
  
  covs <- model.matrix(newformula,data)

  probs <- allProbs(m)
  trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs),m$conditions$betaRef)

  allStates <- NULL
  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(data$ID==unique(data$ID)[zoo])) # nb of observations for animal zoo

    if(zoo!=nbAnimals) {
      p <- probs[aInd[zoo]:(aInd[zoo+1]-1),]
      tm <- trMat[,,aInd[zoo]:(aInd[zoo+1]-1)]
    }
    else {
      p <- probs[aInd[zoo]:nrow(probs),]
      tm <- trMat[,,aInd[zoo]:nrow(probs)]
    }

    xi <- matrix(NA,nbObs,nbStates)
    foo <- (delta[zoo,]%*%tm[,,1])*p[1,]
    xi[1,] <- foo/sum(foo)
    for(i in 2:nbObs) {
      foo <- apply(xi[i-1,]*tm[,,i],2,max)*p[i,]
      xi[i,] <- foo/sum(foo)
    }

    stSeq <- rep(NA,nbObs)
    stSeq[nbObs] <- which.max(xi[nbObs,])
    for(i in (nbObs-1):1)
      stSeq[i] <- which.max(tm[,stSeq[i+1],i+1]*xi[i,])

    allStates <- c(allStates,stSeq)
  }

  return(allStates)
}

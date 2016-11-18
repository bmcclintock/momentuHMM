
#' Viterbi algorithm
#'
#' For a given model, reconstructs the most probable states sequence,
#' using the Viterbi algorithm.
#'
#' @param m An object \code{moveHMM}
#'
#' @return The sequence of most probable states.
#'
#' @examples
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
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

  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  beta <- m$mle$beta
  delta <- m$mle$delta
  stepDist <- m$conditions$stepDist
  angleDist <- m$conditions$angleDist
  omegaDist <- m$conditions$omegaDist
  dryDist <- m$conditions$dryDist
  diveDist <- m$conditions$diveDist
  iceDist <- m$conditions$iceDist
  landDist <- m$conditions$landDist
  stepPar <- m$mle$stepPar
  anglePar <- m$mle$anglePar
  omegaPar <- m$mle$omegaPar
  dryPar <- m$mle$dryPar
  divePar <- m$mle$divePar
  icePar <- m$mle$icePar
  landPar <- m$mle$landPar
  zeroInflation <- m$conditions$zeroInflation

  if(nbStates==1)
    stop("No states to decode (nbStates=1)")

  # identify covariates
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle" & names(data)!="omega" & names(data)!="dry" & names(data)!="dive" & names(data)!="ice" & names(data)!="land")
  nbCovs <- length(covsCol)-1 # substract intercept column
  covs <- data[,covsCol]

  allProbs <- allProbs(data,nbStates,stepDist,angleDist,omegaDist,dryDist,diveDist,iceDist,landDist,stepPar,anglePar,omegaPar,dryPar,divePar,icePar,landPar,zeroInflation)
  trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))

  nbAnimals <- length(unique(data$ID))
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])

  allStates <- NULL
  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(data$ID==unique(data$ID)[zoo])) # nb of observations for animal zoo
    obsInd <- which(!is.na(data$step) & !is.na(data$angle))

    if(zoo!=nbAnimals) {
      p <- allProbs[aInd[zoo]:(aInd[zoo+1]-1),]
      tm <- trMat[,,aInd[zoo]:(aInd[zoo+1]-1)]
    }
    else {
      p <- allProbs[aInd[zoo]:nrow(allProbs),]
      tm <- trMat[,,aInd[zoo]:nrow(allProbs)]
    }

    xi <- matrix(NA,nbObs,nbStates)
    foo <- delta*p[1,]
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


#' Mixture probabilities
#'
#' For a fitted model, this function computes the probability of each individual being in a particular mixture
#'
#' @param m A \code{momentuHMM} object.
#'
#' @return The matrix of individual mixture probabilities, with element [i,j] the probability
#' of individual i being in mixture j
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
#' mixtureProbs(m.mix)
#' }
#' @references
#' Maruotti, A., and T. Ryden. 2009. A semiparametric approach to hidden Markov models under longitudinal observations. Statistics and Computing 19: 381-393.
#'
#' @export
mixtureProbs <- function(m){
  
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")
  
  nbAnimals <- length(unique(m$data$ID))
  mixtures <- m$conditions$mixtures
  
  #if(mixtures==1)
  #  stop("No mixtures to assign probabilities (mixtures=1)")
  
  la <- logAlpha(m) # forward log-probabilities
  
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,max(which(m$data$ID==unique(m$data$ID)[i])))
  
  if(mixtures>1) pie <- m$mle$pi
  else pie <- matrix(1,nbAnimals,1)
  
  # get probability that individual is in each mixture
  mixProbs <- lnum <- matrix(0,nbAnimals,mixtures)
  for(i in 1:nbAnimals){
    pInd <- which(mapply(function(x) isTRUE(all.equal(x,0)),pie))
    if(length(pInd)){
      pie[i,pInd] <- 1.e-100
      pie[i,-pInd] <- pie[i,-pInd] - (1.e-100*length(pInd))/(ncol(pie)-length(pInd))
    }
    for(mix in 1:mixtures){
      c <- max(la[[mix]][aInd[i],]+log(pie[i,mix]))
      lnum[i,mix] <- c + log(sum(exp(la[[mix]][aInd[i],]+log(pie[i,mix])-c)))
    }
    c <- max(lnum[i,])
    mixProbs[i,] <- exp(lnum[i,] - c - log(sum(exp(lnum[i,]-c))))
    mixProbs[i,pInd] <- 0
  }
  colnames(mixProbs) <- paste0("mix",1:mixtures)
  rownames(mixProbs) <- paste0("ID:",unique(m$data$ID))
  
  return(mixProbs)
}
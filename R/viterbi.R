
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

  data <- m$data
  nbStates <- length(m$stateNames)
  beta <- m$mle$beta
  delta <- m$mle$delta
  
  if(nbStates==1)
    stop("No states to decode (nbStates=1)")

  # identify covariates
  formula<-m$conditions$formula
  stateForms<- terms(formula, specials = paste0("state",1:nbStates))
  newformula<-formula
  if(length(unlist(attr(stateForms,"specials")))){
    newForm<-attr(stateForms,"term.labels")[-unlist(attr(stateForms,"specials"))]
    for(i in 1:nbStates){
      if(!is.null(attr(stateForms,"specials")[[paste0("state",i)]])){
        for(j in 1:(nbStates-1)){
          newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("state",i)]]]))
        }
      }
    }
    newformula<-as.formula(paste("~",paste(newForm,collapse="+")))
  }
  formulaStates<-stateFormulas(newformula,nbStates*(nbStates-1),spec="betaCol")
  if(length(unlist(attr(terms(newformula, specials = c(paste0("betaCol",1:(nbStates*(nbStates-1))),"cosinor")),"specials")))){
    allTerms<-unlist(lapply(formulaStates,function(x) attr(terms(x),"term.labels")))
    newformula<-as.formula(paste("~",paste(allTerms,collapse="+")))
    formterms<-attr(terms.formula(newformula),"term.labels")
  } else {
    formterms<-attr(terms.formula(newformula),"term.labels")
    newformula<-formula
  }
  covs <- model.matrix(newformula,data)

  probs <- allProbs(m)
  trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))

  nbAnimals <- length(unique(data$ID))
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])

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
    foo <- (delta%*%tm[,,1])*p[1,]
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

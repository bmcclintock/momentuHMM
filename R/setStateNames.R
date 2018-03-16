
#' Set \code{stateNames} for a \code{momentuHMM}, \code{miHMM}, \code{HMMfits}, or \code{miSum} object

#' @param model \code{\link{momentuHMM}}, \code{miHMM}, \code{\link{HMMfits}}, or \code{\link{miSum}} object
#' @param stateNames Character string providing state names for the model. See \code{\link{fitHMM}} and \code{\link{MIfitHMM}}.
#' 
#' @return \code{model} object with new \code{stateNames} field
#' 
#' @examples
#' m <- example$m
#' mName <- setStateNames(m, stateNames=c("encamped","exploratory"))
#' 
#' @export

setStateNames <- function (model, stateNames) {
  UseMethod("setStateNames")
}

#' @method setStateNames momentuHMM
#' @export
setStateNames.momentuHMM <- function(model, stateNames)
{
  if(!is.character(stateNames) | (length(stateNames)!=length(model$stateNames))) stop("stateNames must be a character vector of length ",length(model$stateNames))
  model$stateNames <- stateNames
  
  # rename columns and rows of MLEs
  for(i in names(model$conditions$dist)){
    if(model$conditions$DMind[[i]]){
      colnames(model$mle[[i]]) <- stateNames
    }
  }
  if(!is.null(model$mle$delta)) colnames(model$mle$delta) <- stateNames
  if(!is.null(model$mle$gamma)){
    colnames(mle$gamma)<-stateNames
    rownames(mle$gamma)<-stateNames
  }
  model$CIreal <- CIreal(model)
  model$CIbeta <- CIbeta(model)
  model
}

#' @method setStateNames miHMM
#' @export
setStateNames.miHMM <- function(model, stateNames)
{

  model$miSum <- setStateNames(model$miSum, stateNames = stateNames)
  
  ind <- which(unlist(lapply(model$HMMfits,is.momentuHMM)))
  nbStates <- unique(unlist(lapply(model$HMMfits[ind],function(x) length(x$stateNames))))
  if(length(nbStates)>1) stop("nbStates differs among the HMMfits models")
  if(!is.character(stateNames) | (length(stateNames)!=nbStates)) stop("stateNames must be a character vector of length ",nbStates)
  for(i in ind){
    model$HMMfits[[i]] <- setStateNames(model$HMMfits[[i]], stateNames = stateNames)
  }  
  model
}

#' @method setStateNames HMMfits
#' @export
setStateNames.HMMfits <- function(model, stateNames)
{
  ind <- which(unlist(lapply(model,is.momentuHMM)))
  nbStates <- unique(unlist(lapply(model[ind],function(x) length(x$stateNames))))
  if(length(nbStates)>1) stop("nbStates differs among the models")
  if(!is.character(stateNames) | (length(stateNames)!=nbStates)) stop("stateNames must be a character vector of length ",nbStates)
  for(i in ind){
    model[[i]] <- setStateNames(model[[i]], stateNames = stateNames)
  }
  model
}

#' @method setStateNames miSum
#' @export
setStateNames.miSum <- function(model, stateNames)
{
  oldStateNames <- model$stateNames
  
  if(!is.character(stateNames) | (length(stateNames)!=length(oldStateNames))) stop("stateNames must be a character vector of length ",length(oldStateNames))
  model$stateNames <- stateNames
  
  for(i in c("beta","real")){
    for(j in names(model$Par[[i]])){
      for(k in names(model$Par[[i]][[j]])){
        for(l in 1:length(stateNames)){
          rownames(model$Par[[i]][[j]][[k]]) <- gsub(oldStateNames[l],stateNames[l],rownames(model$Par[[i]][[j]][[k]]))
          colnames(model$Par[[i]][[j]][[k]]) <- gsub(oldStateNames[l],stateNames[l],colnames(model$Par[[i]][[j]][[k]]))
        }
      }
    }
  }
  for(j in names(model$Par[["timeInStates"]])){
    for(l in 1:length(stateNames)){
      names(model$Par[["timeInStates"]][[j]]) <- gsub(oldStateNames[l],stateNames[l],names(model$Par[["timeInStates"]][[j]]))
    }
  }
  for(j in names(model$Par[["stateProbs"]])){
    for(l in 1:length(stateNames)){
      colnames(model$Par[["stateProbs"]][[j]]) <- gsub(oldStateNames[l],stateNames[l],colnames(model$Par[["stateProbs"]][[j]]))
    }
  }
  model
}

#' Set \code{modelName} for a \code{momentuHMM}, \code{miHMM}, \code{HMMfits}, or \code{miSum} object

#' @param model \code{\link{momentuHMM}}, \code{\link{miHMM}}, \code{\link{HMMfits}}, or \code{\link{miSum}} object
#' @param modelName Character string providing a name for the model. See \code{\link{fitHMM}} and \code{\link{MIfitHMM}}.
#' 
#' @return \code{model} object with new \code{modelName} field
#' 
#' @examples
#' m <- example$m
#' mName <- setModelName(m, modelName="example")
#' 
#' @export

setModelName <- function (model, modelName) {
  if(!is.character(modelName) | length(modelName)>1) stop("modelName must be a character string")
  UseMethod("setModelName")
}

#' @method setModelName momentuHMM
#' @export
setModelName.momentuHMM <- function(model, modelName)
{
  model$modelName <- modelName
  model
}

#' @method setModelName miHMM
#' @export
setModelName.miHMM <- function(model, modelName)
{
  model$miSum$modelName <- modelName
  for(i in 1:length(model$HMMfits)){
    if(is.momentuHMM(model$HMMfits[[i]])) model$HMMfits[[i]]$modelName <- modelName
  }
  model
}

#' @method setModelName HMMfits
#' @export
setModelName.HMMfits <- function(model, modelName)
{
  for(i in 1:length(model)){
    if(is.momentuHMM(model[[i]])) model[[i]]$modelName <- modelName
  }
  model
}

#' @method setModelName miSum
#' @export
setModelName.miSum <- setModelName.momentuHMM
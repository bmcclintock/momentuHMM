
#' Summary \code{moveData}
#' @method summary moveData
#'
#' @param object A \code{moveData} object.
#' @param details \code{TRUE} if quantiles of the covariate values should be printed (default),
#' \code{FALSE} otherwise.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # m is a moveData object (as returned by prepData), automatically loaded with the package
#' data <- example$data
#'
#' summary(data)
#'
#' @export
#' @importFrom stats quantile

summary.moveData <- function(object,details=TRUE,...)
{
  data <- object

  # print animals' IDs and numbers of observations
  nbAnimals <- length(unique(data$ID))
  if(nbAnimals==1)
    cat("Movement data for 1 animal:\n",sep="")
  else
    cat("Movement data for ",nbAnimals," tracks:\n",sep="")
  for(zoo in 1:nbAnimals)
    cat(as.character(unique(data$ID)[zoo])," -- ",
        length(which(data$ID==unique(data$ID)[zoo]))," observations\n",sep="")

  # identify covariates
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle" & names(data)!="states")

  # print covariates' names
  if(length(covsCol)>0) {
    cat("\nCovariate(s): ",sep="")
    if(!details) {
      if(length(covsCol)>1) {
        for(i in 1:(length(covsCol)-1))
          cat(names(data)[covsCol[i]],", ",sep="")
      }
      cat(names(data)[covsCol[length(covsCol)]],"\n")
    } else {
      # print names and quantiles of the covariates
      for(i in 1:length(covsCol)) {
        cat("\n",names(data)[covsCol[i]],"\n")
        q <- quantile(data[,covsCol[i]])
        # insert the mean in the quantiles' vector
        q[6] <- q[5]
        q[5] <- q[4]
        q[4] <- mean(data[,covsCol[i]])
        names(q) <- c("Min.","25%","Median","Mean","75%","Max.")
        print(q)
      }
    }
  } else {
    cat("No covariates.\n")
  }
}

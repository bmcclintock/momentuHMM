
#' Summary \code{momentuHMMData}
#' @method summary momentuHMMData
#'
#' @param object A \code{momentuHMMData} object.
#' @param dataNames Names of the variables to summarize. Default is \code{dataNames=c("step","angle")}.
#' @param animals Vector of indices or IDs of animals for which data will be summarized.
#' Default: \code{NULL} ; data for all animals are summarized.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # data is a momentuHMMData object (as returned by prepData), automatically loaded with the package
#' data <- example$m$data
#'
#' summary(data,dataNames=c("step","angle","cov1","cov2"))
#'
#' @export
#' @importFrom stats quantile

summary.momentuHMMData <- function(object,dataNames=c("step","angle"),animals=NULL,...)
{
    data <- object

    nbAnimals <- length(unique(data$ID))
    #####################################
    ## Define animals to be summarized ##
    #####################################
    if(is.null(animals)) # all animals are summarized
        animalsInd <- 1:nbAnimals
    else {
        if(is.character(animals)) { # animals' IDs provided
          animalsInd <- NULL
          for(zoo in 1:length(animals)) {
            if(length(which(unique(data$ID)==animals[zoo]))==0) # ID not found
              stop("Check animals argument.")

            animalsInd <- c(animalsInd,which(unique(data$ID)==animals[zoo]))
          }
        }

        if(is.numeric(animals)) { # animals' indices provided
          if(length(which(animals<1))>0 | length(which(animals>nbAnimals))>0) # index out of bounds
            stop("Check animals argument.")

          animalsInd <- animals
        }
    }

    # print individuals' IDs and numbers of observations
    nbAnimals <- length(animalsInd)
    if(nbAnimals==1)
        cat("Movement data for 1 individual:\n\n",sep="")
    else
        cat("Movement data for ",nbAnimals," individuals:\n\n",sep="")
    for(zoo in animalsInd)
        cat(as.character(unique(data$ID)[zoo])," -- ",
            length(which(data$ID==unique(data$ID)[zoo]))," observations\n",sep="")

    # identify columns to summarize
    if(any(dataNames %in% names(data))){
        dataNames <- dataNames[which(dataNames %in% names(data))]
    } else stop('dataNames not found in data')
    covsCol <- dataNames

    # print data summaries
    cat("\n\nData summaries:\n\n",sep="")
    print(summary(as.data.frame(data[which(data$ID %in% unique(data$ID)[animalsInd]),covsCol]),digits=2))
}

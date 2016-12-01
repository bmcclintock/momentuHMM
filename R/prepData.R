
#' Preprocessing of the tracking data
#'
#' @param Data A dataframe of data streams, including optionally a field \code{ID}
#' (identifiers for the observed individuals), coordinates from which step length ('step') and turning angle ('angle') area calculated, and \code{covNames} that identify covariates.  If coordinates are included. If step length and turning are to be calculated from coordinates, the \code{coordNames} argument must identify the names for the longitunal and latitudinal coordinates (e.g., \code{coordNames=c("x","y")}).
#' @param type \code{'LL'} if longitude/latitude provided (default), \code{'UTM'} if easting/northing.
#' @param coordNames Names of the columns of coordinates in the data frame. Default: \code{c("x","y")}.
#'
#' @return An object \code{moveData}, i.e. a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{step}{The step lengths - in kilometers if longitude/latitude provided, and in the metrics of
#' the data otherwise}
#' \item{angle}{The turning angles (if any) - in radians}
#' \item{x}{Either easting or longitude}
#' \item{y}{Either norting or latitude}
#' \item{...}{Covariates (if any)}
#'
#' @examples
#' coord1 <- c(1,2,3,4,5,6,7,8,9,10)
#' coord2 <- c(1,1,1,2,2,2,1,1,1,2)
#' Data <- data.frame(coord1=coord1,coord2=coord2)
#' d <- prepData(Data,type='UTM',coordNames=c("coord1","coord2"))
#'
#' @export
#' @importFrom sp spDistsN1

prepData <- function(Data, type=c('LL','UTM'),coordNames=NULL,covNames=NULL)
{
  if(!is.data.frame(Data)) stop("Data must be a data frame")
  distnames<-names(Data)
  if(any(c("step","angle") %in% distnames) & !is.null(coordNames)) stop("Data objects cannot be named 'step' or 'angle';\n  these names are reserved for step lengths and turning angles calculated from coordinates")
  if(!is.null(covNames)) distnames<-distnames[-match(distnames,covNames)]

  if(!is.null(Data$ID)){
    ID <- as.character(Data$ID) # homogenization of numeric and string IDs
    distnames <- distnames[-which(distnames %in% "ID")]
  } 
  else
    ID <- rep("Animal1",length(x)) # default ID if none provided

  if(length(which(is.na(ID)))>0)
    stop("Missing IDs")
  
  # check arguments
  type <- match.arg(type)
  if(!is.null(coordNames)){
    if(length(which(coordNames %in% names(Data)))<2)
      stop("Check the column names of your coordinates.")
    x <- Data[,coordNames[1]]
    y <- Data[,coordNames[2]]
    distnames<-c("step","angle",distnames[-which(distnames %in% coordNames)])
  }

  data <- data.frame(ID=character())

  nbAnimals <- length(unique(ID))

  # check that each animal's observations are contiguous
  for(i in 1:nbAnimals) {
    ind <- which(ID==unique(ID)[i])
    if(length(ind)!=length(ind[1]:ind[length(ind)]))
      stop("Each animal's obervations must be contiguous.")
  }

  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(ID==unique(ID)[zoo]))
    # d = data for one individual
    d <- data.frame(ID=rep(unique(ID)[zoo],nbObs))
    for(j in distnames){
      genData <- rep(NA,nbObs)
      i1 <- which(ID==unique(ID)[zoo])[1]
      i2 <- i1+nbObs-1
      for(i in (i1+1):(i2-1)) {
        if(!is.null(coordNames)){
          if(j=="step" & !is.na(x[i-1]) & !is.na(x[i]) & !is.na(y[i-1]) & !is.na(y[i])) {
            # step length
            genData[i-i1] <- spDistsN1(pts = matrix(c(x[i-1],y[i-1]),ncol=2),
                                      pt = c(x[i],y[i]),
                                      longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
          } else if(j=="angle" & !is.na(x[i-1]) & !is.na(x[i]) & !is.na(x[i+1]) & !is.na(y[i-1]) & !is.na(y[i]) & !is.na(y[i+1])) {
            # turning angle
            genData[i-i1+1] <- turnAngle(c(x[i-1],y[i-1]),
                                       c(x[i],y[i]),
                                       c(x[i+1],y[i+1]))
          } else genData[i-i1] <- Data[[j]][i]
        } else genData[i-i1] <- Data[[j]][i]
      }
      if(j=="step" & !is.null(coordNames)) genData[i2-i1] <- spDistsN1(pts = matrix(c(x[i2-1],y[i2-1]),ncol=2),pt = c(x[i2],y[i2]),longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
      else if(j!="angle" & !is.null(coordNames)) genData[i2-i1] <- Data[[j]][i2]
      d[[j]] <- genData
    }
    data <- rbind(data,d)
  }

  # identify covariate columns
  if(!is.null(covNames)){
    covsCol <- which(names(Data) %in% covNames)
    if(length(covsCol)>0) {
      covs <- Data[,covsCol,drop=FALSE]
      colnames(covs) <- names(Data)[covsCol]
  
      # account for missing values of the covariates
      if(length(which(is.na(covs)))>0)
        warning(paste("There are",length(which(is.na(covs))),
                      "missing covariate values.",
                      "Each will be replaced by the closest available value."))
      for(i in 1:length(covsCol)) {
        if(length(which(is.na(covs[,i])))>0) { # if covariate i has missing values
          if(is.na(covs[1,i])) { # if the first value of the covariate is missing
            k <- 1
            while(is.na(covs[k,i])) k <- k+1
            for(j in k:2) covs[j-1,i] <- covs[j,i]
          }
          for(j in 2:nrow(Data))
            if(is.na(covs[j,i])) covs[j,i] <- covs[j-1,i]
        }
      }
    }
  } else covs<-NULL

  if(!is.null(coordNames)) data <- cbind(data,x=x,y=y)
  if(!is.null(covs)) data <- cbind(data,covs)
  return(moveData(data))
}


#' Preprocessing of the tracking data
#'
#' @param trackData A dataframe of the tracking data, including at least coordinates
#' (either longitude/latitude values or cartesian coordinates), and optionnaly a field \code{ID}
#' (identifiers for the observed individuals). Additionnal fields are considered as covariates.
#' Note that, if the names of the coordinates are not "x" and "y", the \code{coordNames} argument
#' should specified.
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
#' trackData <- data.frame(coord1=coord1,coord2=coord2)
#' d <- prepData(trackData,type='UTM',coordNames=c("coord1","coord2"))
#'
#' @export
#' @importFrom sp spDistsN1

prepData <- function(trackData, type=c('LL','UTM'),coordNames=c("x","y"))
{
  # check arguments
  type <- match.arg(type)
  if(length(which(coordNames %in% names(trackData)))<2)
    stop("Check the columns names of your coordinates.")

  x <- trackData[,coordNames[1]]
  y <- trackData[,coordNames[2]]
  omega_t <- trackData[,"omega"]
  dry_t <- trackData[,"dry"]
  dive_t <- trackData[,"dive"]
  ice_t <- trackData[,"ice"]
  land_t <- trackData[,"land"]

  if(!is.null(trackData$ID))
    ID <- as.character(trackData$ID) # homogenization of numeric and string IDs
  else
    ID <- rep("Animal1",length(x)) # default ID if none provided

  if(length(which(is.na(ID)))>0)
    stop("Missing IDs")

  data <- data.frame(ID=character(),
                     step=numeric(),
                     angle=numeric())

  nbAnimals <- length(unique(ID))

  # check that each animal's observations are contiguous
  for(i in 1:nbAnimals) {
    ind <- which(ID==unique(ID)[i])
    if(length(ind)!=length(ind[1]:ind[length(ind)]))
      stop("Each animal's obervations must be contiguous.")
  }

  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(ID==unique(ID)[zoo]))
    step <- rep(NA,nbObs)
    angle <- rep(NA,nbObs)
    omega <- rep(NA,nbObs)
    dry <- rep(NA,nbObs)
    dive <- rep(NA,nbObs)
    ice <- rep(NA,nbObs)
    land <- rep(NA,nbObs)
    i1 <- which(ID==unique(ID)[zoo])[1]
    i2 <- i1+nbObs-1
    for(i in (i1+1):(i2-1)) {
      if(!is.na(x[i-1]) & !is.na(x[i]) & !is.na(y[i-1]) & !is.na(y[i])) {
        # step length
        step[i-i1] <- spDistsN1(pts = matrix(c(x[i-1],y[i-1]),ncol=2),
                                  pt = c(x[i],y[i]),
                                  longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
        omega[i-i1] <- omega_t[i]
        dry[i-i1] <- dry_t[i]
        dive[i-i1] <- dive_t[i]
        ice[i-i1] <- ice_t[i]
        land[i-i1] <- land_t[i]
      }

      if(!is.na(x[i-1]) & !is.na(x[i]) & !is.na(x[i+1]) & !is.na(y[i-1]) & !is.na(y[i]) & !is.na(y[i+1])) {
        # turning angle
        angle[i-i1+1] <- turnAngle(c(x[i-1],y[i-1]),
                                   c(x[i],y[i]),
                                   c(x[i+1],y[i+1]))
      }
    }
    step[i2-i1] <- spDistsN1(pts = matrix(c(x[i2-1],y[i2-1]),ncol=2),pt = c(x[i2],y[i2]),longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
    omega[i2-i1] <- omega_t[i2]
    dry[i2-i1] <- dry_t[i2]
    dive[i2-i1] <- dive_t[i2]
    ice[i2-i1] <- ice_t[i2]
    land[i2-i1] <- land_t[i2]   
    
    # d = data for one individual
    d <- data.frame(ID=rep(unique(ID)[zoo],nbObs),
                    step=step,
                    angle=angle,
                    omega=omega,
                    dry=dry,
                    dive=dive,
                    ice=ice,
                    land=land)
    # append individual data to output
    data <- rbind(data,d)
  }

  # identify covariate columns
  covsCol <- which(names(trackData)!="ID" & names(trackData)!=coordNames[1] & names(trackData)!=coordNames[2] & names(trackData)!="omega" & names(trackData)!="dry" & names(trackData)!="dive" & names(trackData)!="ice" & names(trackData)!="land")
  if(length(covsCol)>0) {
    covs <- data.frame(trackData[,covsCol]) # to prevent error if nbCovs==1
    colnames(covs) <- names(trackData)[covsCol]

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
        for(j in 2:nrow(trackData))
          if(is.na(covs[j,i])) covs[j,i] <- covs[j-1,i]
      }
    }
  }
  else covs <- NULL

  data <- cbind(data,x=x,y=y)
  if(!is.null(covs)) data <- cbind(data,covs)
  return(moveData(data))
}


#' Preprocessing of the data streams, including calculation of step length, turning angle, and covariates from location data
#'
#' @param Data A dataframe of data streams, including optionally a field \code{ID}
#' (identifiers for the observed individuals), coordinates from which step length ('step') 
#' and turning angle ('angle') area calculated, and any covariates (with names matching \code{covNames}). 
#' If step length and turning angle are to be calculated from coordinates, the \code{coordNames} argument 
#' must identify the names for the longitunal and latitudinal coordinates (e.g., \code{coordNames=c("x","y")}).
#' With the exception of \code{ID} and coordinates, all variables in Data are treated as data streams unless identified
#' as covariates in \code{covNames}.
#' @param type \code{'LL'} if longitude/latitude provided (default), \code{'UTM'} if easting/northing.
#' @param coordNames Names of the columns of coordinates in the \code{Data} data frame. Default: \code{c("x","y")}.
#' @param covNames Names of any covariates in \code{Data} dataframe.
#' are treated as data streams.
#' @param spatialCovs List of raster layer(s) for any spatial covariates not included in \code{Data}. Covariates specified by 'spatialCovs' are
#' extracted based on location data for each time step.

#' @return An object \code{momentuHMMData}
#'
#' @examples
#' coord1 <- c(1,2,3,4,5,6,7,8,9,10)
#' coord2 <- c(1,1,1,2,2,2,1,1,1,2)
#' Data <- data.frame(coord1=coord1,coord2=coord2)
#' d <- prepData(Data,type='UTM',coordNames=c("coord1","coord2"))
#'
#' @export
#' @importFrom sp spDistsN1

prepData <- function(Data, type=c('LL','UTM'),coordNames=NULL,covNames=NULL,spatialCovs=NULL)
{
  if(!is.data.frame(Data)) stop("Data must be a data frame")
  distnames<-names(Data)
  if(any(c("step","angle") %in% distnames) & !is.null(coordNames)) stop("Data objects cannot be named 'step' or 'angle';\n  these names are reserved for step lengths and turning angles calculated from coordinates")
  if(!is.null(covNames)){
    covsCol <- which(distnames %in% covNames)
    if(!length(covsCol)) stop("covNames not found in Data")
    else distnames<-distnames[-covsCol]
  }
  if(!is.null(Data$ID)){
    ID <- as.character(Data$ID) # homogenization of numeric and string IDs
    distnames <- distnames[-which(distnames %in% "ID")]
  } 
  else
    ID <- rep("Animal1",length(x)) # default ID if none provided
  
  if(length(which(is.na(ID)))>0)
    stop("Missing IDs")
  
  if(!is.null(spatialCovs)){
    nbSpatialCovs<-length(names(spatialCovs))
    for(j in 1:nbSpatialCovs)
      if(class(spatialCovs[[j]])!="RasterLayer") stop("spatialCovs must be of class 'RasterLayer'")
    spatialcovnames<-names(spatialCovs)
  } else nbSpatialCovs <- 0
  
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
  
  if(!is.null(coordNames)) {
    data <- cbind(data,x=x,y=y)
    if(nbSpatialCovs){
      spCovs<-numeric()
      xy<-data
      sp::coordinates(xy)<-c("x","y")
      for(j in 1:nbSpatialCovs)
        spCovs<-cbind(spCovs,spatialCovs[[j]][raster::cellFromXY(spatialCovs[[j]],xy)])
      colnames(spCovs)<-spatialcovnames
      data<-cbind(data,spCovs)
    }
  }
  if(!is.null(covs)) data <- cbind(data,covs)
  return(momentuHMMData(data))
}

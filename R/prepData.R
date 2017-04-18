
#' Preprocessing of the data streams and covariates
#' 
#' Preprocessing of the data streams, including calculation of step length, turning angle, and covariates from location data to be suitable for
#' analysis using \code{\link{fitHMM}}
#'
#' @param Data A dataframe of data streams, including optionally a field \code{ID}
#' (identifiers for the observed individuals), coordinates from which step length ('step') 
#' and turning angle ('angle') area calculated, and any covariates (with names matching \code{covNames} and/or \code{angleCovs}). 
#' If step length and turning angle are to be calculated from coordinates, the \code{coordNames} argument 
#' must identify the names for the x- (longitunal) and y- (latitudinal) coordinates.
#' With the exception of \code{ID} and \code{coordNames}, all variables in \code{Data} are treated as data streams unless identified
#' as covariates in \code{covNames} and/or \code{angleCovs}.
#' @param type \code{'UTM'} if easting/northing provided (the default), \code{'LL'} if longitude/latitude.
#' @param coordNames Names of the columns of coordinates in the \code{Data} data frame. Default: \code{c("x","y")}. If \code{coordNames=NULL} then step lengths, turning angles, 
#' and location covariates (i.e., those specified by \code{spatialCovs}, \code{centers}, and \code{angleCovs}) are not calculated.
#' @param covNames Character vector indicating the names of any covariates in \code{Data} dataframe. Any variables in \code{Data} (other than \code{ID}) that are not identified in 
#' \code{covNames} and/or \code{angleCovs} are assumed to be data streams (i.e., missing values will not be accounted for).
#' @param spatialCovs List of \code{\link[raster]{Raster-class}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{\link[raster]{setZ}} and \code{Data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'time').
#' @param centers 2-column matrix providing the x-coordinates (column 1) and y-coordinates (column 2) for any activity centers (e.g., potential 
#' centers of attraction or repulsion) from which distance and angle covariates will be calculated based on the location data. If no row names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'center1.dist', 'center1.angle', 'center2.dist', 'center2.angle'); otherwise the covariate names are derived from the row names
#' of \code{centers} as \code{paste0(rep(rownames(centers),each=2),c(".dist",".angle"))}. As with covariates identified in \code{angleCovs}, note that the angle covariates for each activity center are calculated relative to 
#' the previous movement direction (instead of true north); this is to allow mean turning angle to be modelled as a function of these covariates using circular-circular regression in \code{\link{fitHMM}}
#' or \code{\link{MIfitHMM}}.
#' @param angleCovs Character vector indicating the names of any circular-circular regression covariates in \code{Data} that need conversion from bearing (relative to true north) to turning angle (relative to previous movement direction) 
#' using \code{\link{circAngles}}.

#' @return An object \code{\link{momentuHMMData}}, i.e., a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{...}{Data streams (e.g., 'step', 'angle', etc.)}
#' \item{x}{Either easting or longitude (if \code{coordNames} is specified)}
#' \item{y}{Either norting or latitude (if \code{coordNames} is specified)}
#' \item{...}{Covariates (if any)}
#'
#' @examples
#' coord1 <- c(1,2,3,4,5,6,7,8,9,10)
#' coord2 <- c(1,1,1,2,2,2,1,1,1,2)
#' cov1 <- rnorm(10)
#' 
#' Data <- data.frame(coord1=coord1,coord2=coord2,cov1=cov1)
#' d <- prepData(Data,coordNames=c("coord1","coord2"),covNames="cov1")
#' 
#' # include additional data stream named 'omega'
#' omega <- rbeta(10,1,1)
#' Data <- data.frame(coord1=coord1,coord2=coord2,omega=omega,cov1=cov1)
#' d <- prepData(Data,coordNames=c("coord1","coord2"),covNames="cov1")
#' 
#' # include 'forest' example raster layer as covariate
#' Data <- data.frame(coord1=coord1*1000,coord2=coord2*1000)
#' spatialCov <- list(forest=forest)
#' d <- prepData(Data,coordNames=c("coord1","coord2"),spatialCovs=spatialCov)
#' 
#' # include 2 activity centers
#' Data <- data.frame(coord1=coord1,coord2=coord2,cov1=cov1)
#' d <- prepData(Data,coordNames=c("coord1","coord2"),covNames="cov1",
#'               centers=matrix(c(0,10,0,10),2,2,dimnames=list(c("c1","c2"),NULL)))
#'               
#' # Include angle covariate (relative to true north) that needs conversion to 
#' # turning angle relative to previous movement direction
#' cov2 <- 2*atan(rnorm(10))
#' Data <- data.frame(coord1=coord1,coord2=coord2,cov1=cov1,cov2=cov2)
#' d <- prepData(Data,coordNames=c("coord1","coord2"),covNames="cov1",
#'               angleCovs="cov2")
#' 
#' @export
#' @importFrom sp spDistsN1
#' @importFrom raster cellFromXY getValues getZ

prepData <- function(Data, type=c('UTM','LL'),coordNames=c("x","y"),covNames=NULL,spatialCovs=NULL,centers=NULL,angleCovs=NULL)
{
  if(!is.data.frame(Data)) stop("Data must be a data frame")
  if(any(dim(Data)==0)) stop("Data is empty")
  distnames<-names(Data)
  if(any(c("step","angle") %in% distnames) & !is.null(coordNames)) stop("Data objects cannot be named 'step' or 'angle';\n  these names are reserved for step lengths and turning angles calculated from coordinates")
  if(!is.null(coordNames)){
    if(length(coordNames)!=2) stop('coordNames must be of length 2')
  }
  if(!is.null(covNames) | !is.null(angleCovs)){
    covNames <- unique(c(covNames,angleCovs))
    covsCol <- which(distnames %in% covNames)
    if(!length(covsCol)) stop("covNames not found in Data")
    else distnames<-distnames[-covsCol]
  }
  if(!is.null(Data$ID)){
    ID <- as.character(Data$ID) # homogenization of numeric and string IDs
    distnames <- distnames[-which(distnames %in% "ID")]
  } 
  else
    ID <- rep("Animal1",nrow(Data)) # default ID if none provided
  
  if(length(which(is.na(ID)))>0)
    stop("Missing IDs")
  
  if(!is.null(spatialCovs)){
    if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
    spatialcovnames<-names(spatialCovs)
    if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
    nbSpatialCovs<-length(spatialcovnames)
    for(j in 1:nbSpatialCovs){
      if(!any(class(spatialCovs[[j]]) %in% c("RasterLayer","RasterBrick","RasterStack"))) stop("spatialCovs$",spatialcovnames[j]," must be of class 'RasterLayer', 'RasterStack', or 'RasterBrick'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs$",spatialcovnames[j])
      if(any(class(spatialCovs[[j]]) %in% c("RasterBrick","RasterStack"))){
        if(is.null(raster::getZ(spatialCovs[[j]]))) stop("spatialCovs$",spatialcovnames[j]," is a raster stack or brick that must have set z values (see ?raster::setZ)")
        else if(!(names(attributes(spatialCovs[[j]])$z) %in% names(Data))) stop("spatialCovs$",spatialcovnames[j]," z value '",names(attributes(spatialCovs[[j]])$z),"' not found in Data")
      }
    }
    if(any(spatialcovnames %in% names(Data))) stop("spatialCovs cannot have same names as data")
    if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  } else nbSpatialCovs <- 0
  
  # check arguments
  type <- match.arg(type)
  if(!is.null(coordNames)){
    if(length(which(coordNames %in% names(Data)))<2)
      stop("coordNames not found in Data")
    if(any(names(Data)[which(!(names(Data) %in% coordNames))] %in% c("x","y"))) stop("non-coordinate Data objects cannot be named 'x' or 'y'")
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
    if(!is.null(coordNames))
      if(length(ind)<3) stop('each individual must have at least 3 observations to calculate step lengths and turning angles')
  }
  
  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(ID==unique(ID)[zoo]))
    # d = data for one individual
    d <- data.frame(ID=rep(unique(ID)[zoo],nbObs))
    for(j in distnames){
      genData <- rep(NA,nbObs)
      i1 <- which(ID==unique(ID)[zoo])[1]
      i2 <- i1+nbObs-1
      if(!(j %in% c("step","angle"))){
        genData <- Data[[j]][i1:i2]
      } else if(!is.null(coordNames)){
        for(i in (i1+1):(i2-1)) {
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
          } else genData[i-i1] <- Data[[j]][i-1]
        }
        if(j=="step" & !is.null(coordNames)) {
          genData[i2-i1] <- spDistsN1(pts = matrix(c(x[i2-1],y[i2-1]),ncol=2),pt = c(x[i2],y[i2]),longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
        } 
      } 
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
    class(data$angle)<-c(class(data$angle), "angle")
    if(nbSpatialCovs){
      spCovs<-numeric()
      xy<-data
      sp::coordinates(xy)<-c("x","y")
      for(j in 1:nbSpatialCovs){
        getCells<-raster::cellFromXY(spatialCovs[[j]],xy)
        if(any(is.na(getCells))) stop("Location data are beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
        fullspCovs <- spatialCovs[[j]][getCells]
        if(inherits(spatialCovs[[j]],"RasterLayer")){
          spCovs<-cbind(spCovs,fullspCovs)
        } else {
          tmpspCovs <- numeric(nrow(data))
          zname <- names(attributes(spatialCovs[[j]])$z)
          zvalues <- raster::getZ(spatialCovs[[j]])
          if(!all(unique(Data[[zname]]) %in% zvalues)) stop("Data$",zname," includes z-values with no matching raster layer in spatialCovs$",spatialcovnames[j])
          for(ii in 1:length(zvalues)){
            tmpspCovs[which(Data[[zname]]==zvalues[ii])] <- fullspCovs[which(Data[[zname]]==zvalues[ii]),ii]
          }
          spCovs<-cbind(spCovs,tmpspCovs)
        }
      }
      colnames(spCovs)<-spatialcovnames
      data<-cbind(data,spCovs)
    }
    if(!is.null(centers)){
      if(dim(centers)[2]!=2) stop("centers must be a matrix consisting of 2 columns (i.e., x- and y-coordinates)")
      centerInd <- which(!apply(centers,1,function(x) any(is.na(x))))
      if(length(centerInd)){
        if(is.null(rownames(centers))) centerNames<-paste0("center",rep(centerInd,each=2),".",rep(c("dist","angle"),length(centerInd)))
        else centerNames <- paste0(rep(rownames(centers),each=2),".",rep(c("dist","angle"),length(centerInd)))
        centerCovs <- data.frame(matrix(NA,nrow=nrow(Data),ncol=length(centerInd)*2,dimnames=list(NULL,centerNames)))
        for(zoo in 1:nbAnimals) {
          nbObs <- length(which(ID==unique(ID)[zoo]))
          i1 <- which(ID==unique(ID)[zoo])[1]
          i2 <- i1+nbObs-1
          for(j in 1:length(centerInd)){
            centerCovs[i1,centerNames[(j-1)*2+1:2]]<-distAngle(c(x[i1],y[i1]),c(x[i1],y[i1]),centers[centerInd[j],])
            for(i in (i1+1):i2) {
              centerCovs[i,centerNames[(j-1)*2+1:2]]<-distAngle(c(x[i-1],y[i-1]),c(x[i],y[i]),centers[centerInd[j],])
            }
          }
        }
        data<-cbind(data,centerCovs)
      }  
    }
    if(!is.null(angleCovs)){
      if(!all(angleCovs %in% names(covs))) stop('angleCovs not found in data')
      for(i in angleCovs){
        covs[[i]]<-circAngles(covs[[i]],data)
      }
    }
  }
  if(!is.null(covs)) data <- cbind(data,covs)
  
  data$ID<-as.factor(data$ID)
  
  return(momentuHMMData(data))
}

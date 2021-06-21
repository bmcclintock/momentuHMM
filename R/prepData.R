
#' Preprocessing of the data streams and covariates
#' 
#' Preprocessing of the data streams, including calculation of step length, turning angle, and covariates from location data to be suitable for
#' analysis using \code{\link{fitHMM}}.
#'
#' @param data Either a data frame of data streams or a \code{\link{crwData}} (or \code{\link{crwHierData}}) object (as returned by \code{\link{crawlWrap}}). If \code{data} is a data frame, it can optionally include a field \code{ID}
#' (identifiers for the observed individuals), coordinates from which step length ('step') 
#' and turning angle ('angle') are calculated, and any covariates (with names matching \code{covNames} and/or \code{angleCovs}). 
#' If step length and turning angle are to be calculated from coordinates, the \code{coordNames} argument 
#' must identify the names for the x- (longitunal) and y- (latitudinal) coordinates, and, for hierarchical data, the \code{coordLevel} argument must identify the level of the hierarchy at which the location data are obtained.
#' With the exception of \code{ID}, \code{coordNames}, and, for hierarchical data, \code{level}, all variables in \code{data} are treated as data streams unless identified
#' as covariates in \code{covNames} and/or \code{angleCovs}.
#' @param ... further arguments passed to or from other methods
#' @export
prepData <- function(data, ...) {
  UseMethod("prepData")
}

#' @rdname prepData
#' @method prepData default
#' @param type \code{'UTM'} if easting/northing provided (the default), \code{'LL'} if longitude/latitude. If \code{type='LL'} then step lengths are calculated in kilometers and turning angles are based on initial bearings (see \code{\link{turnAngle}}).
#' Ignored if \code{data} is a \code{\link{crwData}} object.
#' @param coordNames Names of the columns of coordinates in the \code{data} data frame. Default: \code{c("x","y")}. If \code{coordNames=NULL} then step lengths, turning angles, 
#' and location covariates (i.e., those specified by \code{spatialCovs}, \code{centers}, and \code{angleCovs}) are not calculated. Ignored if \code{data} is a \code{\link{crwData}} object.
#' @param covNames Character vector indicating the names of any covariates in \code{data} dataframe. Any variables in \code{data} (other than \code{ID}) that are not identified in 
#' \code{covNames} and/or \code{angleCovs} are assumed to be data streams (i.e., missing values will not be accounted for).
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'time').
#' @param centers 2-column matrix providing the x-coordinates (column 1) and y-coordinates (column 2) for any activity centers (e.g., potential 
#' centers of attraction or repulsion) from which distance and angle covariates will be calculated based on the location data. If no row names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'center1.dist', 'center1.angle', 'center2.dist', 'center2.angle'); otherwise the covariate names are derived from the row names
#' of \code{centers} as \code{paste0(rep(rownames(centers),each=2),c(".dist",".angle"))}. As with covariates identified in \code{angleCovs}, note that the angle covariates for each activity center are calculated relative to 
#' the previous movement direction (instead of standard direction relative to the x-axis); this is to allow the mean turning angle to be modelled as a function of these covariates using circular-circular regression in \code{\link{fitHMM}}
#' or \code{\link{MIfitHMM}}.
#' @param centroids List where each element is a data frame containing the x-coordinates ('x'), y-coordinates ('y'), and times (with user-specified name, e.g., `time') for centroids (i.e., dynamic activity centers where the coordinates can change over time)
#' from which distance and angle covariates will be calculated based on the location data. If any centroids are specified, then \code{data} must include a column indicating the time of each observation, and this column name must match the corresponding user-specified 
#' name of the time column in \code{centroids} (e.g. `time'). Times can be numeric or POSIXt.  If no list names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'centroid1.dist', 'centroid1.angle', 'centroid2.dist', 'centroid2.angle'); otherwise the covariate names are derived from the list names
#' of \code{centroids} as \code{paste0(rep(names(centroids),each=2),c(".dist",".angle"))}. As with covariates identified in \code{angleCovs}, note that the angle covariates for each centroid are calculated relative to 
#' the previous movement direction (instead of standard direction relative to the x-axis); this is to allow the mean turning angle to be modelled as a function of these covariates using circular-circular regression in \code{\link{fitHMM}}
#' or \code{\link{MIfitHMM}}.
#' @param angleCovs Character vector indicating the names of any circular-circular regression angular covariates in \code{data} or \code{spatialCovs} that need conversion from standard direction (in radians relative to the x-axis) to turning angle (relative to previous movement direction) 
#' using \code{\link{circAngles}}.
#' @param altCoordNames Character string indicating an alternative name for the returned location data. If provided, then \code{prepData} will return easting (or longitude) coordinate names as \code{paste0(altCoordNames,".x")} and northing (or latitude) as \code{paste0(altCoordNames,".y")} instead of \code{x} and \code{y}, respectively. This can be useful for location data that are intended to be modeled using a bivariate normal distribution (see \code{\link{fitHMM}}). Ignored unless \code{coordNames} are provided.

#' @return An object \code{\link{momentuHMMData}} or \code{\link{momentuHierHMMData}}, i.e., a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{...}{Data streams (e.g., 'step', 'angle', etc.)}
#' \item{x}{Either easting or longitude (if \code{coordNames} is specified or \code{data} is a \code{crwData} object)}
#' \item{y}{Either norting or latitude (if \code{coordNames} is specified or \code{data} is a \code{crwData} object)}
#' \item{...}{Covariates (if any)}
#' 
#' @details 
#' \itemize{
#' \item If \code{data} is a \code{\link{crwData}} (or \code{\link{crwHierData}}) object, the \code{\link{momentuHMMData}} (or \code{\link{momentuHierHMMData}}) object created by \code{prepData} includes step lengths and turning angles calculated from the best predicted 
#' locations (i.e., \code{crwData$crwPredict$mu.x} and \code{crwData$crwPredict$mu.y}). Prior to using \code{prepData}, additional data streams or covariates unrelated to location (including z-values associated with
#' \code{spatialCovs} raster stacks or bricks) can be merged with the \code{crwData} (or \code{crwHierData}) object using \code{\link{crawlMerge}}.
#' 
#' \item For hierarchical data, \code{data} must include a 'level' field indicating the level of the hierarchy for each observation, and, for location data identified by \code{coordNames}, the \code{coordLevel} argument must indicate the level of the hierarchy at which the location data are obtained.
#' }
#' 
#' @seealso \code{\link{crawlMerge}}, \code{\link{crawlWrap}}, \code{\link{crwData}}
#'
#' @examples
#' coord1 <- c(1,2,3,4,5,6,7,8,9,10)
#' coord2 <- c(1,1,1,2,2,2,1,1,1,2)
#' cov1 <- rnorm(10)
#' 
#' data <- data.frame(coord1=coord1,coord2=coord2,cov1=cov1)
#' d <- prepData(data,coordNames=c("coord1","coord2"),covNames="cov1")
#' 
#' # include additional data stream named 'omega'
#' omega <- rbeta(10,1,1)
#' data <- data.frame(coord1=coord1,coord2=coord2,omega=omega,cov1=cov1)
#' d <- prepData(data,coordNames=c("coord1","coord2"),covNames="cov1")
#' 
#' # include 'forest' example raster layer as covariate
#' data <- data.frame(coord1=coord1*1000,coord2=coord2*1000)
#' spatialCov <- list(forest=forest)
#' d <- prepData(data,coordNames=c("coord1","coord2"),spatialCovs=spatialCov)
#' 
#' # include 2 activity centers
#' data <- data.frame(coord1=coord1,coord2=coord2,cov1=cov1)
#' d <- prepData(data,coordNames=c("coord1","coord2"),covNames="cov1",
#'               centers=matrix(c(0,10,0,10),2,2,dimnames=list(c("c1","c2"),NULL)))
#'               
#' # include centroid
#' data <- data.frame(coord1=coord1,coord2=coord2,cov1=cov1,time=1:10)
#' d <- prepData(data,coordNames=c("coord1","coord2"),covNames="cov1",
#'               centroid=list(centroid=data.frame(x=coord1+rnorm(10),
#'                                                 y=coord2+rnorm(10),
#'                                                 time=1:10)))
#'               
#' # Include angle covariate that needs conversion to 
#' # turning angle relative to previous movement direction
#' u <- rnorm(10) # horizontal component
#' v <- rnorm(10) # vertical component
#' cov2 <- atan2(v,u)
#' data <- data.frame(coord1=coord1,coord2=coord2,cov1=cov1,cov2=cov2)
#' d <- prepData(data,coordNames=c("coord1","coord2"),covNames="cov1",
#'               angleCovs="cov2")
#' 
#' @export
#' @importFrom sp spDistsN1
# @importFrom raster cellFromXY getValues getZ is.factor levels
prepData.default <- function(data, type=c('UTM','LL'), coordNames=c("x","y"), covNames=NULL, spatialCovs=NULL, centers=NULL, centroids=NULL, angleCovs=NULL, altCoordNames=NULL, ...)
{
  if(is.crwData(data)){
    predData <- data$crwPredict
    crwFits <- data$crwFit
    Time.name<-attr(predData,"Time.name")
    if(!is.null(covNames) | !is.null(angleCovs)){
      covNames <- unique(c(covNames,angleCovs[!(angleCovs %in% names(spatialCovs))]))
      if(!length(covNames)) covNames <- NULL
    }
    znames <- unique(unlist(lapply(spatialCovs,function(x) names(attributes(x)$z))))
    if(length(znames))
      if(!all(znames %in% names(predData))) stop("z-values for spatialCovs raster stack or brick not found in ",deparse(substitute(data)),"$crwPredict")
    omitNames <- c("ID","mu.x","mu.y","x","y",covNames,names(spatialCovs),znames)
    if(!is.null(centers)) {
      if(!is.null(rownames(centers))){
        omitNames <- c(omitNames,paste0(rep(rownames(centers),each=2),c(".dist",".angle")))
      } else {
        omitNames <- c(omitNames,paste0("center",rep(1:nrow(centers),each=2),c(".dist",".angle")))
      }
    }
    if(!is.null(centroids)) {
      if(!is.null(names(centroids))){
        omitNames <- c(omitNames,paste0(rep(names(centroids),each=2),c(".dist",".angle")))
      } else {
        omitNames <- c(omitNames,paste0("centroid",rep(1:length(centroids),each=2),c(".dist",".angle")))
      }
    }
    distnames <- names(predData)[which(!(names(predData) %in% omitNames))]
    data <- data.frame(x=predData$mu.x,y=predData$mu.y,predData[,c("ID",distnames,covNames,znames),drop=FALSE])[which(predData$locType=="p" | predData$locType=="f"),]
    type <- 'UTM'
    coordNames <- c("x","y")
  }
  
  if(!is.null(coordNames)){
    if(!is.null(altCoordNames) && (!is.character(altCoordNames) | length(altCoordNames)>1)) stop("'altCoordNames' must be a character string")
    outNames <- c("x","y")
    if(!is.null(altCoordNames)) outNames <- paste0(altCoordNames,".",outNames)
  }
  
  if(!is.data.frame(data)) stop("data must be a data frame")
  if(any(dim(data)==0)) stop("data is empty")
  
  # for directing hierarchical data to prepData.hierarchical
  hierArgs <- list(...)
  if("hierLevels" %in% names(hierArgs)){
    return(prepData.hierarchical(data, type, coordNames, covNames, spatialCovs, centers, centroids, angleCovs, altCoordNames, ...))
  } else if("coordLevel" %in% names(hierArgs)) stop("'coordLevel' cannot be specified unless 'hierLevels' is also specified")
  
  chkDots(...)
  
  distnames<-names(data)
  if(any(c("step","angle") %in% distnames) & !is.null(coordNames)) stop("data objects cannot be named 'step' or 'angle';\n  these names are reserved for step lengths and turning angles calculated from coordinates")
  if(!is.null(coordNames)){
    if(length(coordNames)!=2) stop('coordNames must be of length 2')
  }
  if(!is.null(covNames) | !is.null(angleCovs)){
    covNames <- unique(c(covNames,angleCovs[!(angleCovs %in% names(spatialCovs))]))
    if(length(covNames)){
      covsCol <- which(distnames %in% covNames)
      if(!length(covsCol)) stop("covNames or angleCovs not found in data")
      else distnames<-distnames[-covsCol]
    }
  }
  if(!is.null(data$ID)){
    ID <- as.character(data$ID) # homogenization of numeric and string IDs
    distnames <- distnames[-which(distnames %in% "ID")]
  } 
  else
    ID <- rep("Animal1",nrow(data)) # default ID if none provided
  
  if(length(which(is.na(ID)))>0)
    stop("Missing IDs")
  
  if(!is.null(spatialCovs)){
    if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
    spatialcovnames<-names(spatialCovs)
    if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
    nbSpatialCovs<-length(spatialcovnames)
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package \"raster\" needed for spatial covariates. Please install it.",
           call. = FALSE)
    }
    for(j in 1:nbSpatialCovs){
      if(!inherits(spatialCovs[[j]],c("RasterLayer","RasterBrick","RasterStack"))) stop("spatialCovs$",spatialcovnames[j]," must be of class 'RasterLayer', 'RasterStack', or 'RasterBrick'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs$",spatialcovnames[j])
      if(inherits(spatialCovs[[j]],c("RasterBrick","RasterStack"))){
        if(is.null(raster::getZ(spatialCovs[[j]]))) stop("spatialCovs$",spatialcovnames[j]," is a raster stack or brick that must have set z values (see ?raster::setZ)")
        else if(!(names(attributes(spatialCovs[[j]])$z) %in% names(data))) stop("spatialCovs$",spatialcovnames[j]," z value '",names(attributes(spatialCovs[[j]])$z),"' not found in data")
      }
    }
    if(any(spatialcovnames %in% names(data))) stop("spatialCovs cannot have same names as data")
    if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  } else nbSpatialCovs <- 0
  
  ids <- unique(ID)
  nbAnimals <- length(ids)
  
  # check arguments
  type <- match.arg(type)
  if(type=="LL"){
    if (!requireNamespace("geosphere", quietly = TRUE)) {
      stop("Package \"geosphere\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }
  if(!is.null(coordNames)){
    if(length(which(coordNames %in% names(data)))<2)
      stop("coordNames not found in data")
    if(any(names(data)[which(!(names(data) %in% coordNames))] %in% c("x","y"))) stop("non-coordinate data objects cannot be named 'x' or 'y'")
    x <- data[,coordNames[1]]
    y <- data[,coordNames[2]]
    distnames<-c("step","angle",distnames[-which(distnames %in% coordNames)])
  }
  
  dataHMM <- data.frame(ID=character())
  
  # check that each animal's observations are contiguous
  for(i in 1:nbAnimals) {
    ind <- which(ID==ids[i])
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
      if(!(j %in% c("step","angle")) | is.null(coordNames)){
        genData <- data[[j]][i1:i2]
      } else if(!is.null(coordNames)){
        for(i in (i1+1):(i2-1)) {
          if(j=="step" & !is.na(x[i-1]) & !is.na(x[i]) & !is.na(y[i-1]) & !is.na(y[i])) {
            # step length
            genData[i-i1] <- sp::spDistsN1(pts = matrix(c(x[i-1],y[i-1]),ncol=2),
                                       pt = c(x[i],y[i]),
                                       longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
          } else if(j=="angle" & !is.na(x[i-1]) & !is.na(x[i]) & !is.na(x[i+1]) & !is.na(y[i-1]) & !is.na(y[i]) & !is.na(y[i+1])) {
            # turning angle
            genData[i-i1+1] <- turnAngle(c(x[i-1],y[i-1]),
                                         c(x[i],y[i]),
                                         c(x[i+1],y[i+1]),type)
          }
        }
        if(j=="step" & !is.na(x[i2-1]) & !is.na(x[i2]) & !is.na(y[i2-1]) & !is.na(y[i2])) {
          genData[i2-i1] <- sp::spDistsN1(pts = matrix(c(x[i2-1],y[i2-1]),ncol=2),pt = c(x[i2],y[i2]),longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
        } 
      } 
      d[[j]] <- genData
    }
    dataHMM <- rbind(dataHMM,d)
  }
  
  # identify covariate columns
  if(length(covNames)){
    if(length(covsCol)>0) {
      covs <- data[,covsCol,drop=FALSE]
      colnames(covs) <- names(data)[covsCol]
      
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
          for(j in 2:nrow(data))
            if(is.na(covs[j,i])) covs[j,i] <- covs[j-1,i]
        }
      }
    }
  } else covs<-NULL
  
  if(!is.null(coordNames)) {
    dataHMM[[outNames[1]]] <- x
    dataHMM[[outNames[2]]] <- y
    class(dataHMM$angle)<-c(class(dataHMM$angle), "angle")
    if(nbSpatialCovs){
      spCovs<-numeric()
      xy<-as.data.frame(dataHMM[,outNames])
      sp::coordinates(xy) <- outNames
      for(j in 1:nbSpatialCovs){
        getCells<-raster::cellFromXY(spatialCovs[[j]],xy)
        if(any(is.na(getCells))) stop("Location data are beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
        fullspCovs <- spatialCovs[[j]][getCells]
        if(inherits(spatialCovs[[j]],"RasterLayer")){
          spCovs<-cbind(spCovs,fullspCovs)
        } else {
          tmpspCovs <- numeric(nrow(dataHMM))
          zname <- names(attributes(spatialCovs[[j]])$z)
          if(!(zname %in% names(data))) stop(zname," z-value for ",spatialcovnames[j], "not found in data")
          zvalues <- raster::getZ(spatialCovs[[j]])
          if(!all(unique(data[[zname]]) %in% zvalues)) stop("data$",zname," includes z-values with no matching raster layer in spatialCovs$",spatialcovnames[j])
          for(ii in 1:length(zvalues)){
            tmpspCovs[which(data[[zname]]==zvalues[ii])] <- fullspCovs[which(data[[zname]]==zvalues[ii]),ii]
          }
          spCovs<-cbind(spCovs,tmpspCovs)
        }
      }
      spCovs<-data.frame(spCovs)
      names(spCovs)<-spatialcovnames
      for(j in spatialcovnames){
        if(any(raster::is.factor(spatialCovs[[j]]))){
          spCovs[[j]] <- factor(spCovs[[j]],levels=unique(unlist(raster::levels(spatialCovs[[j]]))))
        }
      }
      #data<-cbind(data,spCovs)
      if(is.null(covs)) covs <- spCovs
      else covs<-cbind(covs,spCovs)
    }
    if(!is.null(centers)){
      if(!is.matrix(centers)) stop("centers must be a matrix")
      if(dim(centers)[2]!=2) stop("centers must be a matrix consisting of 2 columns (i.e., x- and y-coordinates)")
      centerInd <- which(!apply(centers,1,function(x) any(is.na(x))))
      if(length(centerInd)){
        if(is.null(rownames(centers))) centerNames<-paste0("center",rep(centerInd,each=2),".",rep(c("dist","angle"),length(centerInd)))
        else centerNames <- paste0(rep(rownames(centers),each=2),".",rep(c("dist","angle"),length(centerInd)))
        if(any(centerNames %in% names(data))) stop("centers cannot have same names as data")
        centerCovs <- data.frame(matrix(NA,nrow=nrow(data),ncol=length(centerInd)*2,dimnames=list(NULL,centerNames)))
        for(zoo in 1:nbAnimals) {
          nbObs <- length(which(ID==unique(ID)[zoo]))
          i1 <- which(ID==unique(ID)[zoo])[1]
          i2 <- i1+nbObs-1
          for(j in 1:length(centerInd)){
            centerCovs[i1,centerNames[(j-1)*2+1:2]]<-distAngle(c(x[i1],y[i1]),c(x[i1],y[i1]),centers[centerInd[j],],type)
            for(i in (i1+1):i2) {
              centerCovs[i,centerNames[(j-1)*2+1:2]]<-distAngle(c(x[i-1],y[i-1]),c(x[i],y[i]),centers[centerInd[j],],type)
            }
          }
        }
        dataHMM<-cbind(dataHMM,centerCovs)
      }  
    }
    if(!is.null(centroids)){
      if(!is.list(centroids)) stop("centroids must be a named list")
      centroidNames <- character()
      for(j in 1:length(centroids)){
        if(!is.data.frame(centroids[[j]])) stop("each element of centroids must be a data frame")
        if(dim(centroids[[j]])[2]!=3) stop("each element of centroids must be a data frame consisting of 3 columns (i.e., x-coordinate, y-coordinate, and time)")
        if(!all(c("x","y") %in% colnames(centroids[[j]]))) stop("centroids must include 'x' (x-coordinate) and 'y' (y-coordinate) columns")
        if(any(is.na(centroids[[j]]))) stop("centroids cannot contain missing values")
        timeName <- colnames(centroids[[j]])[which(!(colnames(centroids[[j]]) %in% c("x","y")))]
        if(!(timeName %in% names(data))) stop("data must include '",timeName,"' column")
        if(!all(data[,timeName] %in% centroids[[j]][,timeName])) stop("centroids ",timeName," does not span data; each observation time must have a corresponding entry in centroids")
        if(is.null(names(centroids[j]))) centroidNames <- c(centroidNames,paste0("centroid",rep(j,each=2),".",c("dist","angle")))
        else centroidNames <- c(centroidNames,paste0(rep(names(centroids[j]),each=2),".",c("dist","angle")))
        if(any(centroidNames %in% names(data))) stop("centroids cannot have same names as data")
      }
      centroidCovs <- data.frame(matrix(NA,nrow=nrow(data),ncol=length(centroidNames),dimnames=list(NULL,centroidNames)))
      centroidInd <- length(centroidNames)/2
      for(zoo in 1:nbAnimals) {
        nbObs <- length(which(ID==unique(ID)[zoo]))
        i1 <- which(ID==unique(ID)[zoo])[1]
        i2 <- i1+nbObs-1
        for(j in 1:centroidInd){
          centroidCovs[i1,centroidNames[(j-1)*2+1:2]]<-distAngle(c(x[i1],y[i1]),c(x[i1],y[i1]),as.numeric(centroids[[j]][match(data[i1,timeName],centroids[[j]][,timeName]),c("x","y")]),type)
          for(i in (i1+1):i2) {
            centroidCovs[i,centroidNames[(j-1)*2+1:2]]<-distAngle(c(x[i-1],y[i-1]),c(x[i],y[i]),as.numeric(centroids[[j]][match(data[i,timeName],centroids[[j]][,timeName]),c("x","y")]),type)
          }
        }
      }
      dataHMM<-cbind(dataHMM,centroidCovs)
    }
    if(!is.null(angleCovs)){
      if(!all(angleCovs %in% colnames(covs))) stop('angleCovs not found in data or spatialCovs')
      for(i in angleCovs){
        covs[[i]]<-circAngles(covs[[i]],dataHMM,coordNames=outNames)
      }
    }
  }
  if(!is.null(covs)) dataHMM <- cbind(dataHMM,covs)
  
  dataHMM$ID<-as.factor(dataHMM$ID)
  
  if(!is.null(coordNames)) attr(dataHMM,"coords") <- outNames
  
  return(momentuHMMData(dataHMM))
}

#' @rdname prepData
#' @method prepData hierarchical
#' @param hierLevels Character vector indicating the levels of the hierarchy and their order, from top (coarsest scale) to bottom (finest scale), that are included in \code{data$level}. For example, for a 2-level hierarchy then 
#' \code{hierLevels=c("1","2i","2")} indicates \code{data$level} for each observation can be one of three factor levels: "1" (coarse scale), "2i" (initial fine scale), and "2" (fine scale).  Ignored if \code{data} is a \code{\link{crwHierData}} object.
#' @param coordLevel Character string indicating the level of the hierarchy for the location data. If specified, then \code{data} must include a 'level' field indicating the level of the hierarchy for each observation.  Ignored if \code{coordNames} is \code{NULL} or \code{data} is a \code{\link{crwHierData}} object.
#' 
#' @seealso \code{\link{crwHierData}}
#'
#' @export
#' @importFrom sp spDistsN1
# @importFrom raster cellFromXY getValues getZ
prepData.hierarchical <- function(data, type=c('UTM','LL'), coordNames=c("x","y"), covNames=NULL, spatialCovs=NULL, centers=NULL, centroids=NULL, angleCovs=NULL, altCoordNames = NULL, hierLevels, coordLevel, ...)
{
  if(is.crwHierData(data)){
    predData <- data$crwPredict
    crwFits <- data$crwFit
    if(!is.null(covNames) | !is.null(angleCovs)){
      covNames <- unique(c(covNames,angleCovs[!(angleCovs %in% names(spatialCovs))]))
      if(!length(covNames)) covNames <- NULL
    }
    znames <- unique(unlist(lapply(spatialCovs,function(x) names(attributes(x)$z))))
    if(length(znames))
      if(!all(znames %in% names(predData))) stop("z-values for spatialCovs raster stack or brick not found in ",deparse(substitute(data)),"$crwPredict")
    omitNames <- c("ID","mu.x","mu.y","x","y",covNames,names(spatialCovs),znames)
    if(!is.null(centers)) {
      if(!is.null(rownames(centers))){
        omitNames <- c(omitNames,paste0(rep(rownames(centers),each=2),c(".dist",".angle")))
      } else {
        omitNames <- c(omitNames,paste0("center",rep(1:nrow(centers),each=2),c(".dist",".angle")))
      }
    }
    if(!is.null(centroids)) {
      if(!is.null(names(centroids))){
        omitNames <- c(omitNames,paste0(rep(names(centroids),each=2),c(".dist",".angle")))
      } else {
        omitNames <- c(omitNames,paste0("centroid",rep(1:length(centroids),each=2),c(".dist",".angle")))
      }
    }
    distnames <- names(predData)[which(!(names(predData) %in% omitNames))]
    type <- 'UTM'
    coordNames <- c("x","y")
    hierLevels <- levels(data$crwPredict$level)
    coordLevel <- attr(data$crwPredict,"coordLevel")
    data <- data.frame(x=predData$mu.x,y=predData$mu.y,predData[,c("ID",distnames,covNames,znames),drop=FALSE])[which(is.na(predData$locType) | predData$locType!="o"),]
  }
  
  if(!is.null(coordNames)){
    if(!is.null(altCoordNames) && (!is.character(altCoordNames) | length(altCoordNames)>1)) stop("'altCoordNames' must be a character string")
    outNames <- c("x","y")
    if(!is.null(altCoordNames)) outNames <- paste0(altCoordNames,".",outNames)
  }
  
  chkDots(...)
  
  installDataTree()
  
  if(!is.data.frame(data)) stop("data must be a data frame")
  if(any(dim(data)==0)) stop("data is empty")
  distnames<-names(data)#[which(!(names(data) %in% "level"))]
  
  if(is.null(data$level)) stop("data must include a 'level' field")
  if(any(is.na(data$level))) stop("data$level cannot include missing values")
  if(!is.factor(data$level)) data$level <- factor(data$level,levels=hierLevels)
  if(any(is.na(as.numeric(gsub("i","",data$level))))) stop(stop("data$level can only include:\n     '1', '2i', '2', ..., 'Mi', 'M' where M is the number of levels in the hierarchy"))
  if(any(is.na(as.numeric(gsub("i","",hierLevels))))) stop("hierLevels must be ordered factors of the form:\n     '1', '2i', '2', ..., 'Mi', 'M' where M is the number of levels in the hierarchy")
  else {
    maxlevel <- max(as.numeric(gsub("i","",data$level)))
    if(length(hierLevels)!=(2*maxlevel-1)) stop("hierLevels must be of length 2*M-1 where M is the number of levels in the hierarchy")
    checkHierLevels <- character(2*maxlevel-1)
    checkHierLevels[1] <- "1"
    checkHierLevels[seq(2,2*maxlevel-1,2)] <- paste0(seq(2,2*maxlevel-1,2),"i")
    checkHierLevels[seq(3,2*maxlevel-1,2)] <- as.character(seq(2,2*maxlevel-1,2))
    if(!all(hierLevels==checkHierLevels)) stop("hierLevels must be ordered factors of the form:\n     '1', '2i', '2', ..., 'Mi', 'M' where M is the number of levels in the hierarchy")
  }
  if(nlevels(data$level)!=length(hierLevels) || !all(levels(data$level)==hierLevels) || any(is.na(data$level))) stop("data$level not consistent with 'hierLevels'")
  nbLevels <- length(hierLevels)
  for(k in 1:min(nbLevels,nrow(data))){
    if(data$level[k]!=hierLevels[k]) stop("data$level and/or 'hierLevels' not ordered correctly; observation ",k," is level ",data$level[k]," but factor level is ",hierLevels[k])
  }
  
  if(any(c("step","angle") %in% distnames) & !is.null(coordNames)) stop("data objects cannot be named 'step' or 'angle';\n  these names are reserved for step lengths and turning angles calculated from coordinates")
  if(!is.null(coordNames)){
    if(length(coordNames)!=2) stop('coordNames must be of length 2')
    if(!is.character(coordLevel) | length(coordLevel)!=1) stop("coordLevel must be a character string")
    if(!(coordLevel %in% hierLevels)) stop("'coordLevel' must be one of '",paste(hierLevels,collapse="', '"),"'")
  }
  
  
  if(!is.null(covNames) | !is.null(angleCovs)){
    covNames <- unique(c(covNames,angleCovs[!(angleCovs %in% names(spatialCovs))]))
    if(length(covNames)){
      covsCol <- which(distnames %in% covNames)
      if(!length(covsCol)) stop("covNames or angleCovs not found in data")
      else distnames<-distnames[-covsCol]
    }
  }
  if(!is.null(data$ID)){
    ID <- as.character(data$ID) # homogenization of numeric and string IDs
    distnames <- distnames[-which(distnames %in% "ID")]
  } 
  else {
    ID <- rep("Animal1",nrow(data)) # default ID if none provided
    data$ID <- ID
  }
  
  # check order of data$level for each individual
  for(i in unique(data$ID)){
    if(data[which(data$ID==i)[1],"level"]!=hierLevels[1] | data[tail(which(data$ID==i),1),"level"]!=tail(hierLevels,1)) stop("data$level invalid for individual ",i,": level for first observation must be '",hierLevels[1],"'"," and level for last observation must be '",tail(hierLevels,1),"'")
    if(!all(rle(as.character(data[which(data$ID==i),"level"]))$values==hierLevels)) stop("data$level invalid for individual ",i,": pattern not consistent with 'hierLevels'")
  }
  
  if(!is.null(coordNames)) levelData <- data[which(data$level==coordLevel),]
  
  if(length(which(is.na(ID)))>0)
    stop("Missing IDs")
  
  ids <- unique(ID)
  nbAnimals <- length(ids)
  
  if(!is.null(spatialCovs)){
    if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
    spatialcovnames<-names(spatialCovs)
    if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
    nbSpatialCovs<-length(spatialcovnames)
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package \"raster\" needed for spatial covariates. Please install it.",
           call. = FALSE)
    }
    for(j in 1:nbSpatialCovs){
      if(!inherits(spatialCovs[[j]],c("RasterLayer","RasterBrick","RasterStack"))) stop("spatialCovs$",spatialcovnames[j]," must be of class 'RasterLayer', 'RasterStack', or 'RasterBrick'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs$",spatialcovnames[j])
      if(inherits(spatialCovs[[j]],c("RasterBrick","RasterStack"))){
        if(is.null(raster::getZ(spatialCovs[[j]]))) stop("spatialCovs$",spatialcovnames[j]," is a raster stack or brick that must have set z values (see ?raster::setZ)")
        else if(!(names(attributes(spatialCovs[[j]])$z) %in% names(data))) stop("spatialCovs$",spatialcovnames[j]," z value '",names(attributes(spatialCovs[[j]])$z),"' not found in data")
      }
    }
    if(any(spatialcovnames %in% names(data))) stop("spatialCovs cannot have same names as data")
    if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  } else nbSpatialCovs <- 0
  
  # check arguments
  type <- match.arg(type)
  if(type=="LL"){
    if (!requireNamespace("geosphere", quietly = TRUE)) {
      stop("Package \"geosphere\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }
  if(!is.null(coordNames)){
    if(length(which(coordNames %in% names(data)))<2)
      stop("coordNames not found in data")
    if(any(names(data)[which(!(names(data) %in% coordNames))] %in% coordNames)) stop("non-coordinate data objects cannot be named ",coordNames[1]," or ",coordNames[2])
    x <- data[,coordNames[1]]
    y <- data[,coordNames[2]]
    lx <- levelData[,coordNames[1]]
    ly <- levelData[,coordNames[2]]
    dataHMM <- data[,c("ID",distnames[which(!(distnames %in% coordNames))])]
    dataHMM$step <- dataHMM$angle <- rep(NA,nrow(dataHMM))
    distnames<-c("step","angle",distnames[which(!(distnames %in% coordNames))])
  } else {
    dataHMM <- data[,c("ID",distnames)]
  }
  
  # check that each animal's observations are contiguous
  for(i in 1:nbAnimals) {
    ind <- which(data$ID==unique(data$ID)[i])
    if(length(ind)!=length(ind[1]:ind[length(ind)]))
      stop("Each animal's obervations must be contiguous.")
  }
  
  if(!is.null(coordNames)){
    for(zoo in 1:nbAnimals) {
      ind <- which(levelData$ID==unique(levelData$ID)[zoo])
      nbObs <- length(ind)
      if(nbObs<3) stop('each individual must have at least 3 observations to calculate step lengths and turning angles')
      
      # d = data for one individual
      d <- data.frame(ID=rep(unique(levelData$ID)[zoo],nbObs))
      for(j in c("step","angle")){
        genData <- rep(NA,nbObs)
        i1 <- which(levelData$ID==unique(levelData$ID)[zoo])[1]
        i2 <- i1+nbObs-1
        
        for(i in (i1+1):(i2-1)) {
          if(j=="step" & !is.na(lx[i-1]) & !is.na(lx[i]) & !is.na(ly[i-1]) & !is.na(ly[i])) {
            # step length
            genData[i-i1] <- sp::spDistsN1(pts = matrix(c(lx[i-1],ly[i-1]),ncol=2),
                                       pt = c(lx[i],ly[i]),
                                       longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
          } else if(j=="angle" & !is.na(lx[i-1]) & !is.na(lx[i]) & !is.na(lx[i+1]) & !is.na(ly[i-1]) & !is.na(ly[i]) & !is.na(ly[i+1])) {
            # turning angle
            genData[i-i1+1] <- turnAngle(c(lx[i-1],ly[i-1]),
                                         c(lx[i],ly[i]),
                                         c(lx[i+1],ly[i+1]),type)
          }
        }
        if(j=="step" & !is.na(lx[i2-1]) & !is.na(lx[i2]) & !is.na(ly[i2-1]) & !is.na(ly[i2])) {
          genData[i2-i1] <- sp::spDistsN1(pts = matrix(c(lx[i2-1],ly[i2-1]),ncol=2),pt = c(lx[i2],ly[i2]),longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
        }
        dataHMM[which(ID==unique(ID)[zoo] & dataHMM$level==coordLevel),j] <- genData
      } 
    }
  }
  
  # identify covariate columns
  if(length(covNames)){
    if(length(covsCol)>0) {
      covs <- data[,covsCol,drop=FALSE]
      colnames(covs) <- names(data)[covsCol]
      
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
          for(j in 2:nrow(data))
            if(is.na(covs[j,i])) covs[j,i] <- covs[j-1,i]
        }
      }
    }
  } else covs<-NULL
  
  if(!is.null(coordNames)) {
    dataHMM[[outNames[1]]] <- x
    dataHMM[[outNames[2]]] <- y
    if(nbSpatialCovs | !is.null(centers) | !is.null(centroids) | !is.null(angleCovs)){
      # temporarily fill in location data for other levels of the hierarchy to get spatial covariates
      for(zoo in 1:nbAnimals){
        iInd <- which(ID==unique(ID)[zoo])
        X <- dataHMM[[outNames[1]]][iInd]
        X[1] <- X[which.max(!is.na(X))]
        Y <- dataHMM[[outNames[2]]][iInd]
        Y[1] <- Y[which.max(!is.na(Y))]
        for(k in 1:(length(iInd)-1)){
          if(is.na(X[k+1])){
            X[k+1] <- X[k]
          }
          if(is.na(Y[k+1])){
            Y[k+1] <- Y[k]
          }
        }
        dataHMM[[outNames[1]]][iInd] <- X
        dataHMM[[outNames[2]]][iInd] <- Y
      }
      x <- dataHMM[[outNames[1]]]
      y <- dataHMM[[outNames[2]]]
    }
    class(dataHMM$angle)<-c(class(dataHMM$angle), "angle")
    if(nbSpatialCovs){
      spCovs<-numeric()
      xy<-as.data.frame(dataHMM[,outNames])
      sp::coordinates(xy) <- outNames
      for(j in 1:nbSpatialCovs){
        getCells<-raster::cellFromXY(spatialCovs[[j]],xy)
        if(any(is.na(getCells))) stop("Location data are beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
        fullspCovs <- spatialCovs[[j]][getCells]
        if(inherits(spatialCovs[[j]],"RasterLayer")){
          spCovs<-cbind(spCovs,fullspCovs)
        } else {
          tmpspCovs <- numeric(nrow(dataHMM))
          zname <- names(attributes(spatialCovs[[j]])$z)
          if(!(zname %in% names(data))) stop(zname," z-value for ",spatialcovnames[j], "not found in data")
          zvalues <- raster::getZ(spatialCovs[[j]])
          if(!all(unique(data[[zname]]) %in% zvalues)) stop("data$",zname," includes z-values with no matching raster layer in spatialCovs$",spatialcovnames[j])
          for(ii in 1:length(zvalues)){
            tmpspCovs[which(data[[zname]]==zvalues[ii])] <- fullspCovs[which(data[[zname]]==zvalues[ii]),ii]
          }
          spCovs<-cbind(spCovs,tmpspCovs)
        }
      }
      spCovs<-data.frame(spCovs)
      names(spCovs)<-spatialcovnames
      for(j in spatialcovnames){
        if(any(raster::is.factor(spatialCovs[[j]]))){
          spCovs[[j]] <- factor(spCovs[[j]],levels=unique(unlist(raster::levels(spatialCovs[[j]]))))
        }
      }
      #data<-cbind(data,spCovs)
      if(is.null(covs)) covs <- spCovs
      else covs<-cbind(covs,spCovs)
    }
    if(!is.null(centers)){
      if(!is.matrix(centers)) stop("centers must be a matrix")
      if(dim(centers)[2]!=2) stop("centers must be a matrix consisting of 2 columns (i.e., x- and y-coordinates)")
      centerInd <- which(!apply(centers,1,function(x) any(is.na(x))))
      if(length(centerInd)){
        if(is.null(rownames(centers))) centerNames<-paste0("center",rep(centerInd,each=2),".",rep(c("dist","angle"),length(centerInd)))
        else centerNames <- paste0(rep(rownames(centers),each=2),".",rep(c("dist","angle"),length(centerInd)))
        if(any(centerNames %in% names(data))) stop("centers cannot have same names as data")
        centerCovs <- data.frame(matrix(NA,nrow=nrow(data),ncol=length(centerInd)*2,dimnames=list(NULL,centerNames)))
        for(zoo in 1:nbAnimals) {
          nbObs <- length(which(ID==unique(ID)[zoo]))
          i1 <- which(ID==unique(ID)[zoo])[1]
          i2 <- i1+nbObs-1
          for(j in 1:length(centerInd)){
            centerCovs[i1,centerNames[(j-1)*2+1:2]]<-distAngle(c(x[i1],y[i1]),c(x[i1],y[i1]),centers[centerInd[j],],type)
            for(i in (i1+1):i2) {
              centerCovs[i,centerNames[(j-1)*2+1:2]]<-distAngle(c(x[i-1],y[i-1]),c(x[i],y[i]),centers[centerInd[j],],type)
            }
          }
        }
        dataHMM<-cbind(dataHMM,centerCovs)
      }  
    }
    if(!is.null(centroids)){
      if(!is.list(centroids)) stop("centroids must be a named list")
      centroidNames <- character()
      for(j in 1:length(centroids)){
        if(!is.data.frame(centroids[[j]])) stop("each element of centroids must be a data frame")
        if(dim(centroids[[j]])[2]!=3) stop("each element of centroids must be a data frame consisting of 3 columns (i.e., x-coordinate, y-coordinate, and time)")
        if(!all(c("x","y") %in% colnames(centroids[[j]]))) stop("centroids must include 'x' (x-coordinate) and 'y' (y-coordinate) columns")
        if(any(is.na(centroids[[j]]))) stop("centroids cannot contain missing values")
        timeName <- colnames(centroids[[j]])[which(!(colnames(centroids[[j]]) %in% c("x","y")))]
        if(!(timeName %in% names(data))) stop("data must include '",timeName,"' column")
        if(!all(data[,timeName] %in% centroids[[j]][,timeName])) stop("centroids ",timeName," does not span data; each observation time must have a corresponding entry in centroids")
        if(is.null(names(centroids[j]))) centroidNames <- c(centroidNames,paste0("centroid",rep(j,each=2),".",c("dist","angle")))
        else centroidNames <- c(centroidNames,paste0(rep(names(centroids[j]),each=2),".",c("dist","angle")))
        if(any(centroidNames %in% names(data))) stop("centroids cannot have same names as data")
      }
      centroidCovs <- data.frame(matrix(NA,nrow=nrow(data),ncol=length(centroidNames),dimnames=list(NULL,centroidNames)))
      centroidInd <- length(centroidNames)/2
      for(zoo in 1:nbAnimals) {
        nbObs <- length(which(ID==unique(ID)[zoo]))
        i1 <- which(ID==unique(ID)[zoo])[1]
        i2 <- i1+nbObs-1
        for(j in 1:centroidInd){
          centroidCovs[i1,centroidNames[(j-1)*2+1:2]]<-distAngle(c(x[i1],y[i1]),c(x[i1],y[i1]),as.numeric(centroids[[j]][match(data[i1,timeName],centroids[[j]][,timeName]),c("x","y")]),type)
          for(i in (i1+1):i2) {
            centroidCovs[i,centroidNames[(j-1)*2+1:2]]<-distAngle(c(x[i-1],y[i-1]),c(x[i],y[i]),as.numeric(centroids[[j]][match(data[i,timeName],centroids[[j]][,timeName]),c("x","y")]),type)
          }
        }
      }
      dataHMM<-cbind(dataHMM,centroidCovs)
    }
    if(!is.null(angleCovs)){
      if(!all(angleCovs %in% colnames(covs))) stop('angleCovs not found in data or spatialCovs')
      for(i in angleCovs){
        covs[[i]]<-circAngles(covs[[i]],dataHMM,coordNames=outNames)
      }
    }
  }
  if(!is.null(covs)) dataHMM <- cbind(dataHMM,covs)
  
  dataHMM$ID<-as.factor(dataHMM$ID)
  
  if(!is.null(coordNames)) {
    # return coordNames to NA if not at coordLevel
    dataHMM[which(dataHMM$level!=coordLevel),outNames] <- NA
    
    attr(dataHMM,"coords") <- outNames
    attr(dataHMM,"coordLevel") <- coordLevel
  }
  
  return(momentuHierHMMData(dataHMM))
}



#' Preprocessing of the hierarchical data streams and covariates
#' 
#' Preprocessing of the hierarchical data streams, including calculation of step length, turning angle, and covariates from location data to be suitable for
#' analysis using \code{\link{fitHierHMM}}
#'
#' @param data Either a data frame of data streams or a \code{\link{crwHierData}} object (as returned by \code{\link{crawlWrap}}). If \code{data} is a data frame, it can optionally include a field \code{ID}
#' (identifiers for the observed individuals), coordinates from which step length ('step') 
#' and turning angle ('angle') are calculated, and any covariates (with names matching \code{covNames} and/or \code{angleCovs}). 
#' If step length and turning angle are to be calculated from coordinates, the \code{coordNames} argument 
#' must identify the names for the x- (longitunal) and y- (latitudinal) coordinates.
#' With the exception of \code{ID} and \code{coordNames}, all variables in \code{data} are treated as data streams unless identified
#' as covariates in \code{covNames} and/or \code{angleCovs}.
#' @param type \code{'UTM'} if easting/northing provided (the default), \code{'LL'} if longitude/latitude. If \code{type='LL'} then step lengths are calculated in kilometers and turning angles are based on initial bearings (see \code{\link{turnAngle}}).
#' Ignored if \code{data} is a \code{\link{crwHierData}} object.
#' @param coordNames Names of the columns of coordinates in the \code{data} data frame. Default: \code{c("x","y")}. If \code{coordNames=NULL} then step lengths, turning angles, 
#' and location covariates (i.e., those specified by \code{spatialCovs}, \code{centers}, and \code{angleCovs}) are not calculated. Ignored if \code{data} is a \code{\link{crwHierData}} object.
#' @param coordLevel Character string indicating the level of the hierarchy for the location data. Ignored if \code{data} is a \code{\link{crwHierData}} object.
#' @param covNames Character vector indicating the names of any covariates in \code{data} dataframe. Any variables in \code{data} (other than \code{ID}) that are not identified in 
#' \code{covNames} and/or \code{angleCovs} are assumed to be data streams (i.e., missing values will not be accounted for).
#' @param spatialCovs List of \code{\link[raster]{Raster-class}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{\link[raster]{setZ}} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'time').
#' @param centers 2-column matrix providing the x-coordinates (column 1) and y-coordinates (column 2) for any activity centers (e.g., potential 
#' centers of attraction or repulsion) from which distance and angle covariates will be calculated based on the location data. If no row names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'center1.dist', 'center1.angle', 'center2.dist', 'center2.angle'); otherwise the covariate names are derived from the row names
#' of \code{centers} as \code{paste0(rep(rownames(centers),each=2),c(".dist",".angle"))}. As with covariates identified in \code{angleCovs}, note that the angle covariates for each activity center are calculated relative to 
#' the previous movement direction (instead of standard direction relative to the x-axis); this is to allow the mean turning angle to be modelled as a function of these covariates using circular-circular regression in \code{\link{fitHierHMM}}
#' or \code{\link{MIfitHierHMM}}.
#' @param centroids List where each element is a data frame containing the x-coordinates ('x'), y-coordinates ('y'), and times (with user-specified name, e.g., `time') for centroids (i.e., dynamic activity centers where the coordinates can change over time)
#' from which distance and angle covariates will be calculated based on the location data. If any centroids are specified, then \code{data} must include a column indicating the time of each observation, and this column name must match the corresponding user-specified 
#' name of the time column in \code{centroids} (e.g. `time'). Times can be numeric or POSIXt.  If no list names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'centroid1.dist', 'centroid1.angle', 'centroid2.dist', 'centroid2.angle'); otherwise the covariate names are derived from the list names
#' of \code{centroids} as \code{paste0(rep(names(centroids),each=2),c(".dist",".angle"))}. As with covariates identified in \code{angleCovs}, note that the angle covariates for each centroid are calculated relative to 
#' the previous movement direction (instead of standard direction relative to the x-axis); this is to allow the mean turning angle to be modelled as a function of these covariates using circular-circular regression in \code{\link{fitHierHMM}}
#' or \code{\link{MIfitHierHMM}}.
#' @param angleCovs Character vector indicating the names of any circular-circular regression angular covariates in \code{data} or \code{spatialCovs} that need conversion from standard direction (in radians relative to the x-axis) to turning angle (relative to previous movement direction) 
#' using \code{\link{circAngles}}.

#' @return An object \code{\link{momentuHierHMMData}}, i.e., a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{level}{The level of the hierarchy for each observation}
#' \item{...}{Data streams (e.g., 'step', 'angle', etc.)}
#' \item{x}{Either easting or longitude (if \code{coordNames} is specified or \code{data} is a \code{crwHierData} object)}
#' \item{y}{Either norting or latitude (if \code{coordNames} is specified or \code{data} is a \code{crwHierData} object)}
#' \item{...}{Covariates (if any)}
#' 
#' If \code{data} is a \code{\link{crwHierData}} object, the \code{\link{momentuHierHMMData}} object created by \code{prepHierData} includes step lengths and turning angles calculated from the best predicted 
#' locations (\code{crwHierData$crwHierPredict$mu.x} and \code{crwHierData$crwHierPredict$mu.y}). Prior to using \code{prepHierData}, additional data streams or covariates unrelated to location (including z-values associated with
#' \code{spatialCovs} raster stacks or bricks) can be merged with the \code{crwHierData} object using \code{\link{crawlMerge}}.
#' 
#' @seealso \code{\link{crawlMerge}}, \code{\link{crawlWrap}}, \code{\link{crwHierData}}
#'
#' @export
#' @importFrom sp spDistsN1
#' @importFrom raster cellFromXY getValues getZ

prepHierData <- function(data, type=c('UTM','LL'),coordNames=c("x","y"),coordLevel,covNames=NULL,spatialCovs=NULL,centers=NULL,centroids=NULL,angleCovs=NULL)
{
  if(is.crwHierData(data)){
    predData <- data$crwHierPredict
    if(!is.null(covNames) | !is.null(angleCovs)){
      covNames <- unique(c(covNames,angleCovs[!(angleCovs %in% names(spatialCovs))]))
      if(!length(covNames)) covNames <- NULL
    }
    znames <- unique(unlist(lapply(spatialCovs,function(x) names(attributes(x)$z))))
    if(length(znames))
      if(!all(znames %in% names(predData))) stop("z-values for spatialCovs raster stack or brick not found in ",deparse(substitute(data)),"$crwHierPredict")
    omitNames <- c("ID","x","y",covNames,names(spatialCovs),znames)
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
    coordLevel <- attr(data$crwHierPredict,"coordLevel")
    data <- data.frame(x=predData$mu.x,y=predData$mu.y,predData[,c("ID",distnames,covNames,znames),drop=FALSE])[which(is.na(predData$locType) | predData$locType!="o"),]
  }
  if(!is.data.frame(data)) stop("data must be a data frame")
  if(any(dim(data)==0)) stop("data is empty")
  distnames<-names(data)#[which(!(names(data) %in% "level"))]
  
  if(is.null(data$level)) stop("data must include a 'level' field")
  if(!is.factor(data$level)) stop("'level' field must be a factor")
  if(!is.character(coordLevel) | length(coordLevel)!=1) stop("coordLevel must be a character string")
  if(!(coordLevel %in% levels(data$level))) stop("'coordLevel' not found in 'level' field")
  levelData <- data[which(data$level==coordLevel),]
  
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
    for(j in 1:nbSpatialCovs){
      if(!any(class(spatialCovs[[j]]) %in% c("RasterLayer","RasterBrick","RasterStack"))) stop("spatialCovs$",spatialcovnames[j]," must be of class 'RasterLayer', 'RasterStack', or 'RasterBrick'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs$",spatialcovnames[j])
      if(any(class(spatialCovs[[j]]) %in% c("RasterBrick","RasterStack"))){
        if(is.null(raster::getZ(spatialCovs[[j]]))) stop("spatialCovs$",spatialcovnames[j]," is a raster stack or brick that must have set z values (see ?raster::setZ)")
        else if(!(names(attributes(spatialCovs[[j]])$z) %in% names(data))) stop("spatialCovs$",spatialcovnames[j]," z value '",names(attributes(spatialCovs[[j]])$z),"' not found in data")
      }
    }
    if(any(spatialcovnames %in% names(data))) stop("spatialCovs cannot have same names as data")
    if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  } else nbSpatialCovs <- 0
  
  # check arguments
  type <- match.arg(type)
  if(!is.null(coordNames)){
    if(length(which(coordNames %in% names(data)))<2)
      stop("coordNames not found in data")
    if(any(names(data)[which(!(names(data) %in% coordNames))] %in% c("x","y"))) stop("non-coordinate data objects cannot be named 'x' or 'y'")
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

  nbAnimals <- length(unique(ID))
  
  # check that each animal's observations are contiguous
  for(i in 1:nbAnimals) {
    ind <- which(levelData$ID==unique(levelData$ID)[i])
    if(length(ind)!=length(ind[1]:ind[length(ind)]))
      stop("Each animal's obervations must be contiguous.")
    if(!is.null(coordNames))
      if(length(ind)<3) stop('each individual must have at least 3 observations to calculate step lengths and turning angles')
  }
  
  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(levelData$ID==unique(levelData$ID)[zoo]))
    # d = data for one individual
    d <- data.frame(ID=rep(unique(levelData$ID)[zoo],nbObs))
    for(j in c("step","angle")){
      genData <- rep(NA,nbObs)
      i1 <- which(levelData$ID==unique(levelData$ID)[zoo])[1]
      i2 <- i1+nbObs-1
      if(!is.null(coordNames)){
        for(i in (i1+1):(i2-1)) {
          if(j=="step" & !is.na(lx[i-1]) & !is.na(lx[i]) & !is.na(ly[i-1]) & !is.na(ly[i])) {
            # step length
            genData[i-i1] <- spDistsN1(pts = matrix(c(lx[i-1],ly[i-1]),ncol=2),
                                       pt = c(lx[i],ly[i]),
                                       longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
          } else if(j=="angle" & !is.na(lx[i-1]) & !is.na(lx[i]) & !is.na(lx[i+1]) & !is.na(ly[i-1]) & !is.na(ly[i]) & !is.na(ly[i+1])) {
            # turning angle
            genData[i-i1+1] <- turnAngle(c(lx[i-1],ly[i-1]),
                                         c(lx[i],ly[i]),
                                         c(lx[i+1],ly[i+1]),type)
          }
        }
        if(j=="step" & !is.null(coordNames)) {
          genData[i2-i1] <- spDistsN1(pts = matrix(c(lx[i2-1],ly[i2-1]),ncol=2),pt = c(lx[i2],ly[i2]),longlat = (type=='LL')) # TRUE if 'LL', FALSE otherwise
        } 
      } 
      dataHMM[which(ID==unique(ID)[zoo] & dataHMM$level==coordLevel),j] <- genData
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
    dataHMM <- cbind(dataHMM,x=x,y=y)
    if(nbSpatialCovs | !is.null(centers) | !is.null(centroids) | !is.null(angleCovs)){
      # fill in location data for other levels of the hierarchy
      for(zoo in 1:nbAnimals){
        iInd <- which(ID==unique(ID)[zoo])
        X <- dataHMM$x[iInd]
        X[1] <- X[which.max(!is.na(X))]
        Y <- dataHMM$y[iInd]
        Y[1] <- Y[which.max(!is.na(Y))]
        for(k in 1:(length(iInd)-1)){
          if(is.na(X[k+1])){
            X[k+1] <- X[k]
          }
          if(is.na(Y[k+1])){
            Y[k+1] <- Y[k]
          }
        }
        dataHMM$x[iInd] <- X
        dataHMM$y[iInd] <- Y
      }
      x <- dataHMM$x
      y <- dataHMM$y
    }
    class(dataHMM$angle)<-c(class(dataHMM$angle), "angle")
    if(nbSpatialCovs){
      spCovs<-numeric()
      xy<-dataHMM
      sp::coordinates(xy)<-c("x","y")
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
        covs[[i]]<-circAngles(covs[[i]],dataHMM)
      }
    }
  }
  if(!is.null(covs)) dataHMM <- cbind(dataHMM,covs)
  
  dataHMM$ID<-as.factor(dataHMM$ID)
  
  return(momentuHierHMMData(dataHMM))
}

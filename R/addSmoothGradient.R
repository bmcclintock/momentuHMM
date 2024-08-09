#' Add smoothed gradient terms to data
#' 
#' Adds smoothed gradients to \code{\link{momentuHMMData}} or \code{\link{momentuHierHMMData}} objects, which can help reduce bias in gradient-based habitat selection coefficients
#' 
#' @param data A \code{\link{momentuHMMData}} object (as returned by \code{\link{prepData}} or \code{\link{simData}}).
# #' @param CT Logical indicating whether or not \code{data} is for continuous-time models. Default: \code{FALSE} (discrete time).
# #' @param Time.name Character string indicating name of the time column in \code{data} (for continuous-time HMMs). Default: 'time'. Ignored if \code{CT=FALSE}.
# #' @param Time.unit Character string indicating units for time difference between observations (e.g. 'auto', 'secs', 'mins', 'hours', 'days', 'weeks'). Ignored unless \code{data[[Time.name]]} is of class \code{\link[base]{date-time}} or \code{\link[base]{date}}. Default: 'auto'. 
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'time').
#' @param weights Numeric vector indicating the weight for the points in a smoothed version of the gradient, where \code{length(weights)}\eqn{\in \{5,9\}} defines the number of points and \code{sum(weights)=1}. Default: \code{c(1/2,1/8,1/8,1/8,1/8)}, where the first weight (1/2) corresponds to the gradient at the observed location, with the remaining weights (1/8, 1/8, 1/8, 1/8) corresponding to clockwise points around the observed location (i.e. northeast, southeast, southwest, northwest). The smoothed gradients are returned with ``\code{.xs}'' (easting gradient) and ``\code{.ys}'' (northing gradient) suffixes added to the names of \code{spatialCovs}.
#' @param scale Numeric scalar indicating the distance from the observed location to each smoothing point (before adjusting by time step for continuous-time HMMs). Default: \code{NULL}, in which case the square root of the mean squared step lengths (divided by the time step for continuous-time HMMs) across all individuals is used, that is, \code{scale}\eqn{=\sqrt{0.5*mean(step^2/dt)}}, where \eqn{step} is step length and \eqn{dt} is the time step (note \eqn{dt=1} for discrete-time HMMs). The smoothing points are placed at a distance of \code{scale}\eqn{\sqrt{2 dt}} from the observed locations.

#' @return A \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} object, with the smooth gradients added for each spatial covariate.
#' 
#' @references
#' 
#' Blackwell, P. G. and J. Matthiopoulos. 2024. Joint inference for telemetry and spatial survey data. Ecology
#' 
#' @export
#' @importFrom utils head
addSmoothGradient <- function(data,#CT=FALSE,Time.name="time",Time.unit='auto',
                              spatialCovs,weights=c(1/2,1/8,1/8,1/8,1/8),scale=NULL){
  if(!is.momentuHMMData(data)) stop("data must be a momentuHMMData object (as returned by prepData or simData) ")
  if(inherits(data,"ctds")) stop("data cannot be a 'ctds' object (as returned by prepCTDS")
  
  if(isTRUE(attr(data,"CT"))){
    CT <- TRUE
    Time.name <- attr(data,"Time.name")
    Time.unit <- attr(data,"Time.unit")
  } else CT <- FALSE
  
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
    if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  } else stop("spatialCovs must be provided")
  
  if(!is.numeric(weights)) stop("weights must be numeric")
  if(!length(weights) %in% c(5,9)) stop("The length of 'weights' must be 5 or 9")
  if(sum(weights)!=1) stop("weights must sum to 1")
  if(any(paste0(rep(c(spatialcovnames),each=2),c(".xs",".ys")) %in% names(data))) stop(paste0(paste0(rep(c(spatialcovnames),each=2),c(".xs",".ys"))[which(paste0(rep(c(spatialcovnames),each=2),c(".xs",".ys")) %in% names(data))],collapse=", ")," cannot already be in data")
  if(!all(paste0(rep(c(spatialcovnames),each=2),c(".x",".y")) %in% names(data))) stop(paste0(paste0(rep(c(spatialcovnames),each=2),c(".x",".y"))[which(!paste0(rep(c(spatialcovnames),each=2),c(".x",".y")) %in% names(data))],collapse=", ")," gradients not found in data. Was the 'gradient' argument in prepData set to TRUE?")
  
  coordNames <- attr(data,"coords")
  if(is.null(coordNames)) stop("data is missing 'coords' attribute")
  
  npoints <- length(weights)
  
  nbAnimals <- length(unique(data$ID))
  
  # First animal
  iInd <- which(data$ID==unique(data$ID)[1])
  rx <- diff(data[iInd,coordNames[1]])
  ry <- diff(data[iInd,coordNames[2]])
  if(isTRUE(CT)){
    if(is.null(Time.name)) stop("'Time.name' cannot be NULL when CT=TRUE")
    else {
      if(!Time.name %in% names(data)) stop("Time.name not found in 'data'") 
      if("dt" %in% names(data)) stop("'dt' is reserved and cannot be used as a field name in 'data'; please choose a different name") 
      if(inherits(data[[Time.name]],"POSIXt")) {
        tmpTime <- difftime(data[iInd[-1],Time.name],data[iInd[-length(iInd)],Time.name],units=Time.unit)
        Time.unit <- units(tmpTime)
        dt <- as.numeric(tmpTime)
      } else {
        dt <- diff(data[iInd,Time.name])
      }
    }
  } else {
    dt <- rep(1,length(rx))
  }
  locs <- utils::head(data[iInd,coordNames],-1)
  
  # Other animals
  if(nbAnimals>1)
    for (k in 2:nbAnimals)
    {
      iInd <- which(data$ID==unique(data$ID)[k])
      rx <- c(rx,diff(data[iInd,coordNames[1]]))
      ry <- c(ry,diff(data[iInd,coordNames[2]]))
      if(!is.null(Time.name)){
        if(inherits(data[[Time.name]],"POSIXt")){
          dt <- c(dt,as.numeric(difftime(data[iInd[-1],Time.name],data[iInd[-length(iInd)],Time.name],units=Time.unit)))
        } else {
          dt <- c(dt,diff(data[iInd,Time.name]))
        }
      } else dt <- rep(1,length(rx))
      locs <- rbind(locs,head(data[iInd,coordNames],-1))
    }
  
  ddt <- cbind(dt,dt)
  nbObs <- dim(locs)[1]
  
  message("Extracting spatial covariates and calculating smoothed gradients...")
  
  if(is.null(scale)) {
    scale <- sqrt(0.5*mean((rx^2+ry^2)/dt,na.rm=TRUE))
    if(!is.numeric(scale) || length(scale)>1 || scale<0 || !is.finite(scale)) stop("valid scale parameter could not be calculated")
    message("\n    => scale parameter = ",round(scale,6)," \n")
  }
  else {
    if(!is.numeric(scale) || length(scale)>1 || scale<0 || !is.finite(scale)) stop("'scale' must be a positive scalar")
  }
    
  distance <- scale * sqrt(ddt)
  
  points <- vector('list',npoints-1)
  
  newdata <- data[,c("ID",coordNames,paste0(rep(c(spatialcovnames),each=2),c(".x",".y"))[which(paste0(rep(c(spatialcovnames),each=2),c(".x",".y")) %in% names(data))]),drop=FALSE]
  
  if(npoints==9) { # NNE, NEE, SEE, SSE, SSW, SWW, NWW, NNW
    ord <- c(1,(npoints-1):2)
  } else ord <- (npoints-1):1 # NE, SE, SW, NW
  
  oct <- mapply(function(x){
    R <- matrix(c(cos(pi/(npoints-1)),-sin(pi/(npoints-1)),sin(pi/(npoints-1)),cos(pi/(npoints-1))),2,2,byrow=TRUE)
    P <- t(cbind(cos(2*1:(npoints-1)*pi/(npoints-1)),sin(2*1:(npoints-1)*pi/(npoints-1))))*distance[x,1]*sqrt(2);
    pts <- t(R%*%P+matrix(c(locs[x,1],locs[x,2]),nrow=2,ncol=(npoints-1)));
    return(pts[ord,])#[order(mapply(function(y) atan2(pts[y,1]-locs[x,1],pts[y,2]-locs[x,2]),1:nrow(pts))),])
  } ,1:nrow(locs),SIMPLIFY = FALSE)
  
  for(i in 1:(npoints-1)){
    points[[i]] <- t(mapply(function(x) oct[[x]][i,],1:length(oct))) 
  }
  
  obsInd <- unlist(mapply(function(x) utils::head(which(data$ID==x),-1),unique(data$ID),SIMPLIFY = FALSE))
  collapseRast <- tmpColRast <- list()
  for(i in 1:nbSpatialCovs){
    if(inherits(spatialCovs[[i]],c("RasterStack","RasterBrick"))){
      zname <- names(attributes(spatialCovs[[i]])$z)
      zvalues <- raster::getZ(spatialCovs[[i]])
    }
    data[[paste0(spatialcovnames[i],".xs")]] <- newdata[[paste0(spatialcovnames[i],".x")]]
    data[[paste0(spatialcovnames[i],".ys")]] <- newdata[[paste0(spatialcovnames[i],".y")]]
    data[[paste0(spatialcovnames[i],".xs")]][obsInd] <- weights[1] * newdata[[paste0(spatialcovnames[i],".x")]][obsInd]
    data[[paste0(spatialcovnames[i],".ys")]][obsInd] <- weights[1] * newdata[[paste0(spatialcovnames[i],".y")]][obsInd]
    for(j in 1:(npoints-1)){
      for(k in 1:length(coordNames)){
        newdata[[paste0("z",coordNames[k],"_",j)]] <- NA
        newdata[[paste0("z",coordNames[k],"_",j)]][obsInd] <- points[[j]][,k]
      }
      if(inherits(spatialCovs[[i]],"RasterLayer")){
        collapseRast[[spatialcovnames[i]]] <- collapseRaster(spatialCovs[[spatialcovnames[i]]])
        gradtmp <- getGradients(newdata,collapseRast=collapseRast[spatialcovnames[i]],coordNames=c(paste0("z",coordNames[1],"_",j),paste0("z",coordNames[2],"_",j)))
        data[[paste0(spatialcovnames[i],".xs")]][obsInd] <- data[[paste0(spatialcovnames[i],".xs")]][obsInd] + weights[j+1] * gradtmp[[paste0(spatialcovnames[i],".x")]][obsInd]
        data[[paste0(spatialcovnames[i],".ys")]][obsInd] <- data[[paste0(spatialcovnames[i],".ys")]][obsInd] + weights[j+1] * gradtmp[[paste0(spatialcovnames[i],".y")]][obsInd]
      } else if(inherits(spatialCovs[[i]],c("RasterStack","RasterBrick"))){
        collapseRast[[spatialcovnames[i]]] <- list()
        for(ii in 1:length(zvalues)){
          collapseRast[[spatialcovnames[i]]][[which(zvalues==zvalues[ii])]] <- collapseRaster(spatialCovs[[spatialcovnames[i]]][[which(zvalues==zvalues[ii])]])
        }
        tmpColRast[[spatialcovnames[i]]] <- collapseRast[[spatialcovnames[i]]][[which(zvalues==data[[zname]])]]
        gradtmp <- getGradients(newdata,collapseRast=tmpColRast,coordNames=c(paste0("z",coordNames[1],"_",j),paste0("z",coordNames[2],"_",j)))
        data[[paste0(spatialcovnames[i],".xs")]] <- data[[paste0(spatialcovnames[i],".xs")]] + weights[j+1] * gradtmp[[paste0(spatialcovnames[i],".x")]]
        data[[paste0(spatialcovnames[i],".ys")]] <- data[[paste0(spatialcovnames[i],".ys")]] + weights[j+1] * gradtmp[[paste0(spatialcovnames[i],".y")]]
      } 
    }
    progBar(i,nbSpatialCovs)
  }
  
  if(isTRUE(CT)){
    attr(data,"CT") <- TRUE
    attr(data,"Time.name") <- Time.name
    attr(data,"Time.unit") <- Time.unit
  }
  attr(data,"smoothGradient") <- TRUE
  attr(data,"npoints") <- npoints
  return(data)
}
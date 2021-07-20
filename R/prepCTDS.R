
#' Preprocessing of continuous-time discrete-space (CTDS) movement HMMs using ctmcmove
#' 
#' This wrapper function for \code{\link[ctmcmove]{path2ctmc}} and \code{\link[ctmcmove]{ctmc2glm}} converts a \code{data.frame} of coordinates, other data streams, and non-spatial covariates to a \code{\link{momentuHMMData}} object that can be passed directly to \code{\link{fitCTHMM}} (or as a list to \code{\link{MIfitCTHMM}}).
#' 
#' @param data Either a \code{data.frame} of data streams or a \code{\link{crwData}} (or \code{\link{crwHierData}}) object (as returned by \code{\link{crawlWrap}}). If a \code{data.frame}, it must include entries for the x coordinate (\code{x}), y coordinate (\code{y}), and time stamp (\code{time}). An \code{ID} entry must also be included if \code{data} includes multiple individuals.
#' @param ... further arguments passed to or from other methods
#' @export
prepCTDS <- function(data, ...) {
  UseMethod("prepCTDS")
}

#' @rdname prepCTDS
#' @method prepCTDS default
#' @param Time.unit Character string indicating units for time difference between observations (e.g. 'auto', 'secs', 'mins', 'hours', 'days', 'weeks'). Ignored unless \code{data$time} is of class \code{\link[base]{date-time}} or \code{\link[base]{date}}. Default: 'auto', but note that if there are multiple individuals, then the units are determined based on the time stamps for the first individual.
#' @param rast A raster object or raster stack object that will define the discrete-space grid cells for the CTMC movement path.
#' @param directions Integer. Either 4 (indicating a "Rook's neighborhood" of 4 neighboring grid cells) or 8 (indicating a "King's neighborhood" of 8 neighboring grid cells).
#' @param zero.idx Integer vector of the indices of raster cells that are not passable and should be excluded. These are cells where movement should be impossible. Default is zero.idx=integer().
#' @param print.iter Logical. If true, then the progress stepping through each observed location in \code{data} will be output in the terminal.
#' @param interpMethod Specifies interpolation method. Either "ShortestPath", which uses the shortest graphical path on the raster graph, or "LinearInterp", which linearly interpolates between observed locations. "ShortestPath" is slower, slightly more accurate, and allows for impassible barriers specified through "zero.idx". "LinearInterp" is faster but does not allow for impassible barriers.
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'Date'). 
# #' In the \code{\link{momentuHMMData}} object returned by \code{prepCTDS}, covariates for the current position (e.g.\ for use in \code{formula} or \code{DM}) are named with a \code{.cur} suffix (e.g. \code{cov1.cur}).
#' @param spatialCovs.grad List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates, where a directional gradient is to be calculated internally using \code{\link[ctmcmove]{rast.grad}}. Gradient-based covariates specified by \code{spatialCovs.grad} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs.grad} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'Date').
#' @param crw	Logical. If TRUE (default), an autocovariate is created for each cell visited in the CTMC movement path. The autocovariate is a unit-length directional vector pointing from the center of the most recent grid cell to the center of the current grid cell.
#' @param normalize.gradients	Logical. Default is FALSE. If TRUE, then all gradient covariates for \code{spatialCovs.grad} are normalized by dividing by the length of the gradient vector at each point.
#' @param grad.point.decreasing	Logical. Default is FALSE. If TRUE, then the gradient covariates are positive in the direction of decreasing values of the covariate. If FALSE, then the gradient covariates are positive in the direction of increasing values of the covariate (like a true gradient).
#' @param covNames Character vector indicating the names of any covariates in \code{data}. Any variables in \code{data} (other than \code{ID}, \code{x}, \code{y}, and \code{time}) that are not identified in covNames are assumed to be additional data streams (i.e., missing values will not be accounted for).
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing).
#' @return A \code{\link{momentuHMMData}} object of class \code{ctds} that can be passed to \code{\link{fitCTHMM}} (or as a list to \code{\link{MIfitCTHMM}}), where:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{time}{Time stamp for each observation}
#' \item{x}{Easting coordinate of the current cell}
#' \item{y}{Northing coordinate of the current cell}
#' \item{z}{Categorical CTDS data stream indicating which cell was moved to, where \code{z=directions+1} indicates no movement from the current cell. When \code{directions}=4, z=1 indicates a move to the left cell, z=2 right cell, z=3 upward cell, and z=4 downward cell. When \code{directions}=8, z=1 indicates top-left cell, z=2 middle-left cell, z=3 bottom-left cell, z=4 top-right cell, z=5 middle-right cell, z=6 bottom-right cell, z=7 top-middle cell, and z=8 bottom-middle cell.}
#' \item{tau}{Time difference (in \code{Time.unit}) between successive observations}
#' \item{cellCross}{Integer indexing move(s) across multiple cells within a time step (if any)}
#' \item{...}{Additional data streams (if any)}
#' \item{...}{Covariates (if any). Covariates are needed for all neighboring cells and are indexed by \code{1,2,...,directions} using a suffix (e.g. \code{cov.1}, \code{cov.2}, \code{cov.3}, and \code{cov.4} for a covariate named 'cov' with \code{directions=4})}
#' @details 
#' \itemize{
#' \item Any moves to adjacent cells (as defined by \code{directions}; see \code{\link[raster]{adjacent}}) between times t and t+1 are assumed to occur at time t+1 (where \code{tau} = t+1 - t), regardless of the length of the time step (\code{tau}) and where in the current cell or in the neighboring cell that the corresponding (continuous-space) locations reside. 
#' This is done to facilitate the temporal alignment of the discrete-space movement process with other data streams and spatio-temporal covariates. The appropriateness of this simplifying assumption will depend on \code{tau}, the scale of movement, and the raster resolution.
#' 
#' \item If there are move(s) across multiple cells within a time step, data stream(s) other than \code{z} (the categorical CTDS data stream) are set to \code{NA} for these time step(s). 
#' These \code{NA} data stream values must be manually set based on the time spent in each cell if they are to be included in subsequent analysis. 
#' Time-varying spatial covariates for these time step(s) are set based on the z-value for the initial cell and, if different, these must also be set manually based on z-value(s) at the time cell(s) were crossed.
#' All such instances are indicated wherever the 'cellCross' field is > 0.
#' 
#' \item If \code{data} is a \code{\link{crwData}} object, the \code{\link{momentuHMMData}} object created by \code{prepCTDS} is based on the best predicted 
#' locations (i.e., \code{crwData$crwPredict$mu.x} and \code{crwData$crwPredict$mu.y}). Prior to using \code{prepCTDS}, additional data streams or covariates unrelated to location (including z-values associated with
#' \code{spatialCovs} and \code{spatialCovs.grad} raster stacks or bricks) can be merged with the \code{crwData} object using \code{\link{crawlMerge}}.
#' }
#' @references
#' 
#' Hanks E. M., Hooten M. B., and Alldredge M. W. 2015. Continuous-time Discrete-space Models for Animal Movement. The Annals of Applied Statistics 9:145-165
#' 
#' @export
prepCTDS.default <- function(data, Time.unit="auto", rast, directions=4, zero.idx=integer(), print.iter=FALSE, interpMethod="ShortestPath",
                     spatialCovs=NULL, spatialCovs.grad=NULL, crw = TRUE, normalize.gradients = FALSE, grad.point.decreasing = FALSE,
                     covNames=NULL, ncores=1, ...) {
  
  if (!requireNamespace("ctmcmove", quietly = TRUE)) {
    stop("Package \"ctmcmove\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package \"raster\" needed for spatial covariates. Please install it.",
           call. = FALSE)
  }
  
  if(is.crwData(data)){
    predData <- data$crwPredict
    crwFits <- data$crwFit
    Time.name<-attr(predData,"Time.name")
    znames <- unique(c(unlist(lapply(spatialCovs,function(x) names(attributes(x)$z))),unlist(lapply(spatialCovs.grad,function(x) names(attributes(x)$z)))))
    if(length(znames))
      if(!all(znames %in% names(predData))) stop("z-values for spatialCovs or spatialCovs.grad raster stack (or brick) not found in ",deparse(substitute(data)),"$crwPredict")
    omitNames <- c("ID","mu.x","mu.y","x","y","locType",Time.name,covNames,names(spatialCovs),names(spatialCovs.grad),znames)
    distnames <- names(predData)[which(!(names(predData) %in% omitNames))]
    data <- data.frame(x=predData$mu.x,y=predData$mu.y,time=predData[[Time.name]],predData[,c("ID",distnames,covNames,znames),drop=FALSE])[which(predData$locType=="p" | predData$locType=="f"),]
    coordNames <- c("x","y")
  }
  
  if(is.null(data$x) | is.null(data$y) | is.null(data$time))
    stop("data must contain 'x', 'y', and 'time' fields")
  
  if(inherits(data$time,"POSIXt")){
    if (!requireNamespace("lubridate", quietly = TRUE)) {
      stop("Package \"lubridate\" needed for POSIXt times. Please install it.",
           call. = FALSE)
    }
  }
  
  if(!(directions %in% c(4,8))) stop("'directions' must be 4 or 8")
  if(any(names(data) %in% c("x.current","y.current","z",paste0("z.",1:directions),"tau","cellCross"))) stop("'x.current', 'y.current', 'z', 'tau', and 'cellCross' are reserved and cannot be fields in data")
  if(!is.null(covNames)){
    if(any(covNames %in% c(c("ID","time","x","y")))) stop("covNames cannot include 'ID', 'time', 'x', or 'y'")
  }
  
  data <- prepData(data,coordNames=NULL,covNames=covNames)
  dataStreams <- names(data)[which(!names(data) %in% c("ID","time","x","y",covNames))]
  
  # check rasters using prepData
  checkRast(data,spatialCovs,spatialCovs.grad)
  
  if(inherits(rast,"RasterBrick")) rast <- raster::stack(rast)
  
  # if 'auto' time units, use that for first individual
  if(inherits(data$time,"POSIXt") & Time.unit=="auto"){
    iInd <- which(data$ID==unique(data$ID)[1])
    Time.unit <- attr(difftime(data$time[iInd[-1]],data$time[iInd[-length(iInd)]],units="auto"),"units")
  }
  
  chkDots(...)
  
  iDat <- id <- NULL # get rid of no visible binding for global variable warning
  if(ncores>1){
    for(pkg in c("doFuture","future")){
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package \"",pkg,"\" needed for parallel processing to work. Please install it.",
             call. = FALSE)
      }
    }
    oldDoPar <- doFuture::registerDoFuture()
    on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
    future::plan(future::multisession, workers = ncores)
    # hack so that foreach %dorng% can find internal momentuHMM variables without using ::: (forbidden by CRAN)
    progBar <- progBar
    pkgs <- c("momentuHMM")
  } else { 
    doParallel::registerDoParallel(cores=ncores)
    pkgs <- c("momentuHMM")
  }
  #ctdsglm <- list()
  #for(i in unique(data$ID)){
  if(!exists("messInd")) messInd <- TRUE # message indicator exported from MIfitCTHMM foreach call
  ids <- unique(data$ID)
  withCallingHandlers(
    ctdsglm <- foreach(iDat=mapply(function(x) data[which(data$ID==x),],unique(data$ID),SIMPLIFY = FALSE), id=ids, .combine = 'rbind', .export="messInd") %dorng% {
      if(messInd){
        if(ncores==1) message("\rIndividual ",which(ids==id),"... ",sep="")
        else progBar(which(ids==id), length(ids))
      }
      ctds <- path2ctds(xy=as.matrix(iDat[which(!is.na(iDat$x) & !is.na(iDat$y)),c("x","y")]),t=iDat$time[which(!is.na(iDat$x) & !is.na(iDat$y))],rast=rast,directions=directions,zero.idx=zero.idx,print.iter=print.iter,method=interpMethod,Time.unit=Time.unit)
      ctdsglm <- ctds2glm(iDat[which(!is.na(iDat$x) & !is.na(iDat$y)),], ctds,rast = rast, directions=directions, spatialCovs = spatialCovs, spatialCovs.grad=spatialCovs.grad, crw=crw, normalize.gradients = normalize.gradients, grad.point.decreasing = grad.point.decreasing, include.cell.locations = TRUE, zero.idx=zero.idx, covNames = covNames)
      ctdsglm$ID <- id
      if(length(dataStreams) | !is.null(covNames)){
        multiCellMove <- which(ctdsglm$cellCross>0)
        if(length(multiCellMove)){
          warning("There were ",length(unique(ctdsglm$cellCross[multiCellMove]))," move(s) across multiple cells within a time step:\n",
          #"   -- any spatial covariates for ",paste0(dataStreams,collapse=", ")," pertain to the initial cell for these time step(s)\n",
          #"   -- during model fitting, the state(s) for ",paste0(dataStreams,collapse=", ")," are assumed to be the initial state(s) at the start of these time step(s)")
          ifelse(length(dataStreams),paste0("   -- '",paste0(dataStreams,collapse="', '"),"' data stream(s) were set to NA for these time step(s)\n"),""),
          ifelse(length(dataStreams),"   -- data stream values must be manually set based on the time spent in each cell if these are to be included in subsequent analysis\n",""),
          ifelse(!is.null(covNames),"   -- time-varying spatial covariates set based on the z-value for the initial cell and, if different, must be set manually based on z-value(s) at the time cell(s) were crossed\n",""),
          "   -- these instances are indicated wherever the 'cellCross' field is > 0")
          #for(j in unique(ctdsglm$cellCross[multiCellMove])){
            #crInd <- which(ctdsglm$cellCross==j)
            #tmp <- ctdsglm[crInd,][1:directions,]
            #tmp$tau <- sum(ctdsglm[which(ctdsglm$cellCross==j),"tau"])/directions
            #tmp$z <- NA
            #ctdsglm[which(ctdsglm$cellCross==j),dataStreams] <- NA
            #ctdsglm <- rbind(ctdsglm[1:(crInd[1]-1),],tmp,ctdsglm[-(1:(crInd[1]-1)),])
          #}
          ctdsglm[which(ctdsglm$cellCross>0),dataStreams] <- NA
        }
        #ctdsglm[(1:nrow(ctdsglm))[-seq(1,nrow(ctdsglm),directions)],dataStreams] <- NA
      }
      names(ctdsglm)[which(names(ctdsglm)=="t")] <- "time"
      return(ctdsglm)
    },
  warning=muffleRNGwarning)
  if(ncores==1) {
    doParallel::stopImplicitCluster()
    if(messInd) cat("DONE\n")
  }
  else future::plan(future::sequential)
  #}
  #ctdsglm <- do.call(rbind,ctdsglm)
  #rownames(ctdsglm) <- NULL
  
  ctdsglm <- ctdsglm[,which(!colnames(ctdsglm) %in% c("x.adj","y.adj"))]
  ctdsOut <- ctdsglm[seq(1,nrow(ctdsglm),directions),]
  names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x.current","y.current","tau","cellCross",dataStreams,covNames))] <- paste0(names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x.current","y.current","tau","cellCross",dataStreams,covNames))],".1")
  for(j in 2:directions){
    tmp <- ctdsglm[seq(j,nrow(ctdsglm),directions),which(!colnames(ctdsglm) %in% c("ID","step","angle","time","x.current","y.current","tau","cellCross",dataStreams,covNames))]
    names(tmp) <- paste0(names(tmp),".",j)
    ctdsOut <- cbind(ctdsOut,tmp)
  }
  
  # add non-gradient spatial covariates for current position (e.g. for inclusion in TPM)
  #if(!is.null(spatialCovs)) names(spatialCovs) <- paste0(names(spatialCovs),".cur")
  #ctdsOut <- prepData(ctdsOut,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs)
  ctdsOut <- prepData(ctdsOut,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs)
  #ctdsOut <- prepData(ctdsglm,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs)
  zMat <- as.matrix(ctdsOut[,paste0("z.",1:directions)])
  zMat <- cbind(zMat,1-rowSums(zMat))
  ctdsOut$z <- NA
  ctdsOut$z[which(!is.na(ctdsOut$z.1))] <- unlist(apply(zMat,1,which.max))
  ctdsOut[paste0("z.",1:directions)] <- NULL
  ctdsOut <- ctdsOut[,c("ID","time","x","y","z",dataStreams,"tau","cellCross",covNames,names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x","y","step","angle","z",dataStreams,"tau","cellCross",covNames))])]
  if(!any(ctdsOut$cellCross>0)) ctdsOut$cellCross <- NULL
  
  class(ctdsOut) <- unique(append(c("momentuHMMData","ctds"),class(ctdsOut)))
  attr(ctdsOut,"Time.unit") <- Time.unit
  attr(ctdsOut,"directions") <- directions
  attr(ctdsOut,"coords") <- c("x","y")
  attr(ctdsOut,"ctdsData") <- "z"
  attr(ctdsOut,"normalize.gradients") <- normalize.gradients 
  attr(ctdsOut,"grad.point.decreasing") <- grad.point.decreasing
  return(ctdsOut)
}

#' @rdname prepCTDS
#' @method prepCTDS hierarchical
#' @param hierLevels Character vector indicating the levels of the hierarchy and their order, from top (coarsest scale) to bottom (finest scale), that are included in \code{data$level}. For example, for a 2-level hierarchy then 
#' \code{hierLevels=c("1","2i","2")} indicates \code{data$level} for each observation can be one of three factor levels: "1" (coarse scale), "2i" (initial fine scale), and "2" (fine scale).  Ignored if \code{data} is a \code{\link{crwHierData}} object.
#' @param coordLevel Character string indicating the level of the hierarchy for the location data. If specified, then \code{data} must include a 'level' field indicating the level of the hierarchy for each observation.  Ignored if \code{coordNames} is \code{NULL} or \code{data} is a \code{\link{crwHierData}} object.
#' 
#' @export
#' @importFrom sp spDistsN1
# @importFrom raster cellFromXY getValues getZ
prepCTDS.hierarchical <- function(data, Time.unit="auto", rast, directions=4, zero.idx=integer(), print.iter=FALSE, interpMethod="ShortestPath",
                                  spatialCovs=NULL, spatialCovs.grad=NULL, crw = TRUE, normalize.gradients = FALSE, grad.point.decreasing = FALSE,
                                  covNames=NULL, ncores=1, hierLevels, coordLevel, ...) {
  
  if (!requireNamespace("ctmcmove", quietly = TRUE)) {
    stop("Package \"ctmcmove\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package \"raster\" needed for spatial covariates. Please install it.",
         call. = FALSE)
  }
  
  if(is.crwHierData(data)){
    predData <- data$crwPredict
    crwFits <- data$crwFit
    Time.name<-attr(predData,"Time.name")
    znames <- unique(unlist(lapply(spatialCovs,function(x) names(attributes(x)$z))))
    if(length(znames))
      if(!all(znames %in% names(predData))) stop("z-values for spatialCovs raster stack or brick not found in ",deparse(substitute(data)),"$crwPredict")
    omitNames <- c("ID","mu.x","mu.y","x","y","locType",Time.name,covNames,names(spatialCovs),names(spatialCovs.grad),znames)
    distnames <- names(predData)[which(!(names(predData) %in% omitNames))]
    #type <- 'UTM'
    coordNames <- c("x","y")
    hierLevels <- levels(data$crwPredict$level)
    coordLevel <- attr(data$crwPredict,"coordLevel")
    data <- data.frame(x=predData$mu.x,y=predData$mu.y,time=predData[[Time.name]],predData[,c("ID",distnames,covNames,znames),drop=FALSE])[which(is.na(predData$locType) | predData$locType!="o"),]
    for(i in unique(data$ID)){
      iInd <- (data$ID==i)
      cInd <- (data$level==coordLevel)
      for(h in hierLevels[which(!hierLevels %in% coordLevel)]){
        hInd <- (data$level==h)
        data[which(iInd & hInd),c("x","y")] <- data[which(iInd & cInd),][match(data$time[which(iInd & hInd)],data$time[which(iInd & cInd)]),c("x","y")]
      }
      iInd <- which(iInd)
      for(j in iInd[-length(iInd)]){
        if(is.na(data[j+1,"x"])) data[j+1,"x"] <- data[j,"x"]
        if(is.na(data[j+1,"y"])) data[j+1,"y"] <- data[j,"y"]
      }
    }
  }
  
  if(is.null(data$x) | is.null(data$y) | is.null(data$time))
    stop("data must contain 'x', 'y', and 'time' fields")
  
  if(inherits(data$time,"POSIXt")){
    if (!requireNamespace("lubridate", quietly = TRUE)) {
      stop("Package \"lubridate\" needed for POSIXt times. Please install it.",
           call. = FALSE)
    }
  }
  
  if(any(names(data) %in% c("x.current","y.current","z",paste0("z.",1:directions),"tau","cellCross"))) stop("'x.current', 'y.current', 'z', 'tau', and 'cellCross' are reserved and cannot be fields in data")
  if(!is.null(covNames)){
    if(any(covNames %in% c(c("ID","time","x","y")))) stop("covNames cannot include 'ID', 'time', 'x', or 'y'")
  }
  
  data <- prepData(data,coordNames=NULL,covNames=covNames, hierLevels=hierLevels,coordLevel=coordLevel)
  dataStreams <- names(data)[which(!names(data) %in% c("ID","time","x","y",covNames))]
  
  # check rasters using prepData
  checkRast(data[which(data$level==coordLevel),],spatialCovs,spatialCovs.grad)
  
  if(inherits(rast,"RasterBrick")) rast <- raster::stack(rast)
  
  # if 'auto' time units, use that for first individual
  if(inherits(data$time,"POSIXt") & Time.unit=="auto"){
    iInd <- which(data$ID==unique(data$ID)[1])
    Time.unit <- attr(difftime(data$time[iInd[-1]],data$time[iInd[-length(iInd)]],units="auto"),"units")
  }
  
  chkDots(...)
  
  iDat <- id <- NULL # get rid of no visible binding for global variable warning
  if(ncores>1){
    for(pkg in c("doFuture","future")){
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package \"",pkg,"\" needed for parallel processing to work. Please install it.",
             call. = FALSE)
      }
    }
    oldDoPar <- doFuture::registerDoFuture()
    on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
    future::plan(future::multisession, workers = ncores)
    # hack so that foreach %dorng% can find internal momentuHMM variables without using ::: (forbidden by CRAN)
    progBar <- progBar
    path2ctds <- path2ctds
    ctds2glm <- ctds2glm
    pkgs <- c("momentuHMM")
  } else { 
    doParallel::registerDoParallel(cores=ncores)
    pkgs <- c("momentuHMM")
  }
  #ctdsglm <- list()
  #for(i in unique(data$ID)){
  if(!exists("messInd")) messInd <- TRUE # message indicator exported from MIfitCTHMM foreach call
  ids <- unique(data$ID)
  withCallingHandlers(
    ctdsglm <- foreach(iDat=mapply(function(x) data[which(data$ID==x),],unique(data$ID),SIMPLIFY = FALSE), id=ids, .combine = 'rbind', .export="messInd") %dorng% {
      if(messInd){
        if(ncores==1) message("\rIndividual ",which(ids==id),"... ",sep="")
        else progBar(which(ids==id), length(ids))
      }
      ctdsCoord <- path2ctds(xy=as.matrix(iDat[which(!is.na(iDat$x) & !is.na(iDat$y) & iDat$level==coordLevel),c("x","y")]),t=iDat$time[which(!is.na(iDat$x) & !is.na(iDat$y) & iDat$level==coordLevel)],rast=rast,directions=directions,zero.idx=zero.idx,print.iter=print.iter,method=interpMethod,Time.unit=Time.unit)
      ctds <- path2ctds(xy=as.matrix(iDat[which(!is.na(iDat$x) & !is.na(iDat$y)),c("x","y")]),t=iDat$time[which(!is.na(iDat$x) & !is.na(iDat$y))],rast=rast,directions=directions,zero.idx=zero.idx,print.iter=print.iter,method=interpMethod,Time.unit=Time.unit)
      lInd <- which(!is.na(iDat$x) & !is.na(iDat$y) & iDat$level==coordLevel)
      ctds$ec[lInd] <- ctdsCoord$ec
      ctds$rt[lInd] <- ctdsCoord$rt
      ctds$trans.times[lInd] <- ctdsCoord$trans.times
      ctds$cellCross[lInd] <- ctdsCoord$cellCross
      for(iLevel in hierLevels){
        #if(iLevel!=coordLevel){
          lInd <- which(iDat$level[which(!is.na(iDat$x) & !is.na(iDat$y))]==iLevel)
          if(!grepl("i",iLevel)) ctds$rt[lInd] <- c(diff(iDat[which(!is.na(iDat$x) & !is.na(iDat$y)),"time"][lInd]),0)
          else ctds$rt[lInd] <- 1
        #}
      }
      class(ctds) <- append("hierarchical",class(ctds))
      attr(ctds,"coordLevel") <- coordLevel
      ctdsglmCoord <- ctds2glm(iDat[which(!is.na(iDat$x) & !is.na(iDat$y) & iDat$level==coordLevel),], ctdsCoord,rast = rast, directions=directions, spatialCovs = spatialCovs, spatialCovs.grad=spatialCovs.grad, crw=crw, normalize.gradients = normalize.gradients, grad.point.decreasing = grad.point.decreasing, include.cell.locations = TRUE, zero.idx=zero.idx, covNames = covNames)
      ctdsglm <- ctds2glm(iDat[which(!is.na(iDat$x) & !is.na(iDat$y)),], ctds,rast = rast, directions=directions, spatialCovs = spatialCovs, spatialCovs.grad=spatialCovs.grad, crw=crw, normalize.gradients = normalize.gradients, grad.point.decreasing = grad.point.decreasing, include.cell.locations = TRUE, zero.idx=zero.idx, covNames = covNames)
      ctdsglm[which(ctdsglm$level==coordLevel & ctdsglm$t %in% ctdsglmCoord$t),] <- ctdsglmCoord
      ctdsglm[which(ctdsglm$level==coordLevel & (!ctdsglm$t %in% ctdsglmCoord$t)),"z"] <- NA
      ctdsglm$ID <- id
      if(length(dataStreams) | !is.null(covNames)){
        multiCellMove <- which(ctdsglm$cellCross>0)
        if(length(multiCellMove)){
          warning("There were ",length(unique(ctdsglm$cellCross[multiCellMove]))," move(s) across multiple cells within a time step:\n",
                  #"   -- any spatial covariates for ",paste0(dataStreams,collapse=", ")," pertain to the initial cell for these time step(s)\n",
                  #"   -- during model fitting, the state(s) for ",paste0(dataStreams,collapse=", ")," are assumed to be the initial state(s) at the start of these time step(s)")
                  ifelse(length(dataStreams),paste0("   -- '",paste0(dataStreams,collapse="', '"),"' data stream(s) were set to NA for these time step(s)\n"),""),
                  ifelse(length(dataStreams),"   -- data stream values must be manually set based on the time spent in each cell if these are to be included in subsequent analysis\n",""),
                  ifelse(!is.null(covNames),"   -- time-varying spatial covariates set based on the z-value for the initial cell and, if different, must be set manually based on z-value(s) at the time cell(s) were crossed\n",""),
                  "   -- these instances are indicated wherever the 'cellCross' field is > 0")
          #for(j in unique(ctdsglm$cellCross[multiCellMove])){
          #crInd <- which(ctdsglm$cellCross==j)
          #tmp <- ctdsglm[crInd,][1:directions,]
          #tmp$tau <- sum(ctdsglm[which(ctdsglm$cellCross==j),"tau"])/directions
          #tmp$z <- NA
          #ctdsglm[which(ctdsglm$cellCross==j),dataStreams] <- NA
          #ctdsglm <- rbind(ctdsglm[1:(crInd[1]-1),],tmp,ctdsglm[-(1:(crInd[1]-1)),])
          #}
          ctdsglm[which(ctdsglm$cellCross>0),dataStreams[which(dataStreams!="level")]] <- NA
        }
        #ctdsglm[(1:nrow(ctdsglm))[-seq(1,nrow(ctdsglm),directions)],dataStreams] <- NA
      }
      names(ctdsglm)[which(names(ctdsglm)=="t")] <- "time"
      return(ctdsglm)
    },
    warning=muffleRNGwarning)
  if(ncores==1) {
    doParallel::stopImplicitCluster()
    if(messInd) cat("DONE\n")
  }
  else future::plan(future::sequential)
  #}
  #ctdsglm <- do.call(rbind,ctdsglm)
  #rownames(ctdsglm) <- NULL
  
  ctdsglm <- ctdsglm[,which(!colnames(ctdsglm) %in% c("x.adj","y.adj"))]
  ctdsOut <- ctdsglm[seq(1,nrow(ctdsglm),directions),]
  names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x.current","y.current","tau","cellCross",dataStreams,covNames))] <- paste0(names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x.current","y.current","tau","cellCross",dataStreams,covNames))],".1")
  for(j in 2:directions){
    tmp <- ctdsglm[seq(j,nrow(ctdsglm),directions),which(!colnames(ctdsglm) %in% c("ID","step","angle","time","x.current","y.current","tau","cellCross",dataStreams,covNames))]
    names(tmp) <- paste0(names(tmp),".",j)
    ctdsOut <- cbind(ctdsOut,tmp)
  }
  
  # add non-gradient spatial covariates for current position (e.g. for inclusion in TPM)
  #if(!is.null(spatialCovs)) names(spatialCovs) <- paste0(names(spatialCovs),".cur")
  #ctdsOut <- prepData(ctdsOut,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs)
  ctdsOut <- prepData(ctdsOut,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs)
  #ctdsOut <- prepData(ctdsglm,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs)
  zMat <- as.matrix(ctdsOut[,paste0("z.",1:directions)])
  zMat <- cbind(zMat,1-rowSums(zMat))
  ctdsOut$z <- NA
  ctdsOut$z[which(!is.na(ctdsOut$z.1))] <- unlist(apply(zMat,1,which.max))
  ctdsOut[paste0("z.",1:directions)] <- NULL
  ctdsOut <- ctdsOut[,c("ID","time","x","y","z",dataStreams,"tau","cellCross",covNames,names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x","y","step","angle","z",dataStreams,"tau","cellCross",covNames))])]
  if(!any(ctdsOut$cellCross>0)) ctdsOut$cellCross <- NULL
  
  ctdsOut$z[which(ctdsOut$level!=coordLevel)] <- NA
  class(ctdsOut) <- unique(append(c("momentuHierHMMData","ctds"),class(ctdsOut)))
  attr(ctdsOut,"Time.unit") <- Time.unit
  attr(ctdsOut,"directions") <- directions
  attr(ctdsOut,"coords") <- c("x","y")
  attr(ctdsOut,"ctdsData") <- "z"
  attr(ctdsOut,"normalize.gradients") <- normalize.gradients 
  attr(ctdsOut,"grad.point.decreasing") <- grad.point.decreasing
  return(ctdsOut)
}

path2ctds<-function (xy, t, rast, directions = 4, zero.idx = integer(), 
          print.iter = FALSE, method = "ShortestPath", Time.unit = "auto") 
{
  if (class(rast) == "RasterStack") {
    rast = rast[[1]]
  }
  if(!(directions %in% c(4,8))) stop("'directions' must be 4 or 8")
  
  raster::values(rast) <- 1
  raster::values(rast)[zero.idx] <- 0
  trans = gdistance::transition(rast, prod, directions = directions)
  ncell = raster::ncell(rast)
  #if (method == "LinearInterp") {
    A = Matrix::Matrix(0, nrow = ncell, ncol = ncell, sparse = TRUE)
    adj = raster::adjacent(rast, 1:ncell, directions=directions)
    A[adj] <- 1
  #}
    
  if(inherits(t,"POSIXt") & Time.unit=="auto"){
    Time.unit <- attr(difftime(t[-1],t[-length(t)],units="auto"),"units")
  }
      
  path = cbind(xy, as.numeric(t-min(t),units=Time.unit))
  tidx = sort(as.numeric(t), index.return = T)$ix
  path = path[tidx, ]
  T = nrow(path)
  ec.all = raster::cellFromXY(rast, xy)
  ec = ec.all[1]
  current.cell = ec
  rt = integer()
  cellCross <- integer()
  current.cr <- 0
  if (print.iter) {
    cat("Total locations =", T, "\n")
  }
  for (i in 2:(T)) {
    if (print.iter) {
      cat(i, " ")
    }
    if (ec.all[i] == current.cell) {
      rt <- c(rt,as.numeric(t[i] - t[i - 1],units=Time.unit))
      ec = c(ec, ec.all[i])
      cellCross <- c(cellCross,0)
      current.rt <- rt[length(rt)]
    }
    else {
      if (A[current.cell, ec.all[i]] == 1) {
        rt = c(rt, as.numeric(t[i] - t[i - 1],units=Time.unit))
        ec = c(ec, ec.all[i])
        current.cell = ec.all[i]
        cellCross <- c(cellCross,0)
      } else if (method == "ShortestPath") {
        sl = gdistance::shortestPath(trans, as.numeric(xy[i - 1, 
        ]), as.numeric(xy[i, ]), "SpatialLines")
        slc = raster::coordinates(sl)[[1]][[1]]
        sl.cells = raster::cellFromXY(rast, slc)
        t.in.each.cell = as.numeric(t[i] - t[i - 1],units=Time.unit)/length(sl.cells)
        rt = c(rt, rep(t.in.each.cell, length(sl.cells)))
        ec = c(ec, sl.cells)
        current.cell = ec[length(ec)]
        current.cr <- current.cr + 1
        cellCross <- c(cellCross,rep(current.cr,length(sl.cells)))
        current.rt <- t.in.each.cell
      } else if (method == "LinearInterp") {
        xyt.1 = path[i - 1, ]
        xyt.2 = path[i, ]
        d = sqrt(sum((xyt.1[-3] - xyt.2[-3])^2))
        rast.res = raster::res(rast)[1]
        xapprox = stats::approx(c(xyt.1[3], xyt.2[3]), c(xyt.1[1], 
                                                  xyt.2[1]), n = max(100, round(d/rast.res * 
                                                                                  100)))
        yapprox = stats::approx(c(xyt.1[3], xyt.2[3]), c(xyt.1[2], 
                                                  xyt.2[2]), n = max(100, round(d/rast.res * 
                                                                                  100)))
        tapprox = xapprox$x
        xapprox = xapprox$y
        yapprox = yapprox$y
        xycell.approx = raster::cellFromXY(rast, cbind(xapprox, 
                                               yapprox))
        rle.out = rle(xycell.approx)
        ec.approx = rle.out$values
        transition.idx.approx = c(1,cumsum(rle.out$lengths))
        transition.times.approx = tapprox[transition.idx.approx]
        rt.approx = diff(transition.times.approx)#diff(c(t[i-1],transition.times.approx))#
        if(!rt.approx[1]>0){
          #ec <- ec[-length(ec)]
          rt.approx <- rt.approx[-1]
          ec.approx <- ec.approx[-1]
          #ec.approx <- c(ec.approx,ec.approx[length(ec.approx)])
        }
        rt = c(rt, rt.approx)
        ec = c(ec, ec.approx)
        current.cell = ec[length(ec)]
        current.cr <- current.cr + 1
        cellCross <- c(cellCross,rep(current.cr,length(rt.approx)))
        current.rt <- rt.approx[length(rt.approx)]
      } else stop("interpMethod must be 'ShortestPath' or 'LinearInterp'")
    }
  }
  if(inherits(t,"POSIXt")) trans.times <- c(t[1],t[1] + lubridate::duration(cumsum(rt),units=Time.unit))
  else trans.times <- c(t[1],t[1] + cumsum(rt))
  list(ec = ec, rt = c(diff(c(path[1,3],path[1,3] + cumsum(rt))),0), trans.times = trans.times, cellCross=c(cellCross,0), Time.unit=Time.unit)
}

ctds2glm <- function (data, ctmc, rast, spatialCovs=NULL, spatialCovs.grad=NULL, crw = TRUE, normalize.gradients = FALSE, 
          grad.point.decreasing = TRUE, include.cell.locations = TRUE, 
          directions = 4, zero.idx = integer(), covNames) 
{
  if(length(spatialCovs)){
    examplerast <- spatialCovs[[1]]
  } else examplerast <- rast
  
  if(!(directions %in% c(4,8))) stop("'directions' must be 4 or 8")
  
  locs = ctmc$ec
  wait.times = ctmc$rt
  notzero.idx = 1:raster::ncell(examplerast)
  if (length(zero.idx) > 0) {
    notzero.idx = notzero.idx[-zero.idx]
  }
  adj = raster::adjacent(examplerast, locs, pairs = TRUE, sorted = TRUE, 
                         id = TRUE, directions = directions, target = notzero.idx)
  
  adjCount <- c(unname(table(adj[,"id"])))
  if(!all(adjCount==directions)) stop("locations cannot reside in cells along the edge of the raster; please expand the extent of 'rast' accordingly")
  
  adj.cells = adj[, 3]
  rr = rle(adj[, 1])
  time.idx = rep(rr$values, times = rr$lengths)
  start.cells = adj[, 2]
  moves <- which(locs[-length(locs)]!=locs[-1])
  z = rep(0, length(start.cells))
  idx.move = rep(0, length(z))
  diag.move = rep(0, length(locs))
  for (i in moves) {
    idx.t = which(time.idx == i)
    idx.m = which(adj.cells[idx.t] == locs[i + 1])
    z[idx.t[idx.m]] <- 1
    if (length(idx.m) == 0) {
      diag.move[i] = 1
      #z[idx.t] <- NA
      stop("diagonal move detected and cannot be properly accounted for",ifelse(directions==4,"; try expanding 'directions' to 8",""))
    }
  }
  tau = rep(wait.times, times = rr$lengths)
  t = rep(ctmc$trans.times, times = rr$lengths)
  cr <- rep(ctmc$cellCross, times = rr$length)
  
  xy.cell = raster::xyFromCell(examplerast, start.cells)
  xy.adj = raster::xyFromCell(examplerast, adj.cells)
  v.adj = (xy.adj - xy.cell)/sqrt(apply((xy.cell - xy.adj)^2, 
                                        1, sum))
  
  if(inherits(ctmc,"hierarchical")){
    moveData <- merge(data.frame(ID=data$ID[1],time=ctmc$trans.times[which(data$level==attr(ctmc,"coordLevel"))]),data,by=c("time","ID"),all=TRUE)
    moveData <- moveData[order(moveData$time,moveData$level),]
  } else {
    moveData <- merge(data.frame(ID=data$ID[1],time=ctmc$trans.times),data,by=c("time","ID"),all=TRUE)
    moveData <- moveData[order(moveData$time),]
  }
  moveData <- suppressWarnings(prepData(moveData,coordNames=NULL,covNames=covNames))
  moveData <- moveData[which(moveData$time %in% ctmc$trans.times),]
  #if(any(order(moveData$time)!=(1:nrow(moveData)))) moveData <- moveData[order(moveData$time),] # check for any weird time sorting by merge
  dataStreams <- names(moveData)[which(!names(moveData) %in% c("ID","time","x","y",covNames))]
  
  ### extract spatial covariates
  # gradient-based covariates
  p.grad = length(spatialCovs.grad)

  if(p.grad) {
    X.grad = get.grad(moveData,covNames,directions,start.cells,v.adj,spatialCovs.grad,normalize.gradients,grad.point.decreasing)
  }
  
  # non-gradient based covariates
  X.static <- get.static(moveData,covNames,directions,start.cells,spatialCovs)
  
  zInd <- which(rowSums(matrix(z,ncol=directions,byrow=TRUE))==1)
  if(!length(zInd)) stop("no moves were made between cells")
  if(zInd[length(zInd)]*directions<length(z)) zInd <- c(zInd,zInd[length(zInd)]+1)
  moveInd = c(t(matrix(1:length(z),ncol=directions,byrow=TRUE)[zInd,]))
  X.crw <- numeric(length(z))
  idx.move = which(z == 1)
  idx.move = c(idx.move, min(length(z),tail(idx.move,1)+directions))
  v.moves = v.adj[rep(idx.move[1:(length(rr$lengths[zInd]) - 1)], 
                      times = rr$lengths[zInd][-1]), ]
  v.moves = rbind(matrix(0, ncol = 2, nrow = directions), 
                  v.moves)
  X.crw[moveInd] = apply(v.moves * v.adj[moveInd,], 1, sum)
  zInd <- c(0,zInd)
  for(i in 1:(length(zInd)-1)){
    if(zInd[i+1]-zInd[i] > 1) {
      X.crw[(zInd[i]*directions+1):((zInd[i+1]-1)*directions)] <- rep(matrix(X.crw,ncol=directions,byrow=TRUE)[zInd[i+1],],(zInd[i+1]-1)-zInd[i])
    }
  }
  X.crw[((zInd[i])*directions+1):length(X.crw)] <- matrix(X.crw,ncol=directions,byrow=TRUE)[zInd[i+1],]
  
  if (crw == FALSE & p.grad > 0) {
    X = cbind(X.static, X.grad)
  }
  if (crw == TRUE & p.grad > 0) {
    X = cbind(X.static, X.grad, X.crw)
    colnames(X)[ncol(X)] = "crw"
  }
  if (crw == FALSE & p.grad == 0) {
    X = cbind(X.static)
  }
  if (crw == TRUE & p.grad == 0) {
    X = cbind(X.static, X.crw)
    colnames(X)[ncol(X)] = "crw"
  }
  if (include.cell.locations) {
    xys = cbind(xy.cell, xy.adj)
    colnames(xys) = c("x.current", "y.current", "x.adj", 
                      "y.adj")
    X = cbind(X, xys)
  }
  if(length(dataStreams)) X = cbind(X,moveData[,dataStreams,drop=FALSE][rep(seq_len(nrow(moveData)),each=directions),,drop=FALSE])
  if(length(covNames)) X = cbind(X,moveData[,covNames,drop=FALSE][rep(seq_len(nrow(moveData)),each=directions),,drop=FALSE])
  out = data.frame(z = z, X, tau = tau, t = t, cellCross = cr)
  T = nrow(out)
  out = out[-((T - (directions-1)):T), ]
  out
}

checkRast <- function(data,spatialCovs,spatialCovs.grad){
  check1 <- prepData(data.frame(data[1:3,]),spatialCovs=spatialCovs)
  check2 <- prepData(data.frame(data[1:3,]),spatialCovs=spatialCovs.grad)
  if(!is.null(spatialCovs)) if(any(is.na(raster::cellFromXY(spatialCovs[[1]],data[,c("x","y")])))) stop("Location data are beyond the spatial extent of the raster(s). Try expanding the extent of the raster(s).")
  else if(!is.null(spatialCovs.grad)) if(any(is.na(raster::cellFromXY(spatialCovs.grad[[1]],data[,c("x","y")])))) stop("Location data are beyond the spatial extent of the raster(s). Try expanding the extent of the raster(s).")
  if(any(names(spatialCovs) %in% names(spatialCovs.grad))) stop("'spatialCovs' and 'spatialCovs.grad' names must be unique")
}

#insertRow <- function(existingDF, newrow, r) {
#  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
#  existingDF[r,] <- newrow
#  existingDF
#}

get.grad <- function(moveData,covNames,directions,start.cells,v.adj,spatialCovs.grad,normalize.gradients,grad.point.decreasing,gradMat=TRUE)
{
  X.grad <- do.call(ifelse(gradMat,"cbind","list"),lapply(spatialCovs.grad,
                       function(x){
                         if(inherits(x,"RasterLayer")){
                           x = ctmcmove::rast.grad(x)
                           if (normalize.gradients) {
                             lengths = sqrt(x$grad.x^2 + x$grad.y^2)
                             x$grad.x <- x$grad.x/lengths
                             x$grad.y <- x$grad.y/lengths
                           }
                           if(any(!is.finite(x$grad.x)) | any(!is.finite(x$grad.y))){
                             warning("some gradients were not finite and were set to zero")
                             if(any(!is.finite(x$grad.x))) x$grad.x[which(!is.finite(x$grad.x))] <- 0
                             if(any(!is.finite(x$grad.y))) x$grad.y[which(!is.finite(x$grad.y))] <- 0
                           }
                           return(v.adj[, 1] * x$grad.x[start.cells] + v.adj[, 2] * x$grad.y[start.cells])
                         } else {
                           zname <- names(attributes(x)$z)
                           zvalues <- raster::getZ(x)
                           if(inherits(x,"RasterBrick")) x <- raster::setZ(raster::stack(x),zvalues,zname)
                           grad <- ctmcmove::rast.grad(x[[1]])
                           if (normalize.gradients) {
                             lengths = sqrt(grad$rast.grad.x^2 + grad$rast.grad.y^2)
                             grad$rast.grad.x <- grad$rast.grad.x/lengths
                             grad$rast.grad.y <- grad$rast.grad.y/lengths
                           }
                           gradx <- raster::stack(grad$rast.grad.x)
                           grady <- raster::stack(grad$rast.grad.y)
                           for(i in 2:raster::nlayers(x)){
                             grad <- ctmcmove::rast.grad(x[[i]])
                             if (normalize.gradients) {
                               lengths = sqrt(grad$rast.grad.x^2 + grad$rast.grad.y^2)
                               grad$rast.grad.x <- grad$rast.grad.x/lengths
                               grad$rast.grad.y <- grad$rast.grad.y/lengths
                             }
                             gradx <- raster::stack(gradx,grad$rast.grad.x)
                             grady <- raster::stack(grady,grad$rast.grad.y)
                           }
                           names(gradx) <- names(grady) <- names(x)
                           gradx <- raster::setZ(gradx,zvalues,zname)
                           grady <- raster::setZ(grady,zvalues,zname)
                           fullx <- gradx[start.cells]
                           fully <- grady[start.cells]
                           if(gradMat){
                             tmpspCovs.x <- tmpspCovs.y <- numeric(length(start.cells))
                             if((!zname %in% covNames) & zname!="time") stop("z-value name '",zname,"' must be included in covNames")
                             for(ii in 1:length(zvalues)){
                               tmpspCovs.x[which(rep(moveData[[zname]],each=directions)==zvalues[ii])] <- fullx[which(rep(moveData[[zname]],each=directions)==zvalues[ii]),ii]
                               tmpspCovs.y[which(rep(moveData[[zname]],each=directions)==zvalues[ii])] <- fully[which(rep(moveData[[zname]],each=directions)==zvalues[ii]),ii]
                             }
                             return(v.adj[, 1] * tmpspCovs.x + v.adj[, 2] * tmpspCovs.y)
                           } else {
                             return(v.adj[, 1] * fullx + v.adj[, 2] * fully)
                           }
                         }
                       }))
  
  if (grad.point.decreasing == TRUE) {
    if(is.list(X.grad)) X.grad <- lapply(X.grad,function(x) -x)
    else X.grad = -X.grad
  }
  return(X.grad)
}

get.static <- function(moveData,covNames,directions,start.cells,spatialCovs){
  do.call(cbind,lapply(spatialCovs,
                       function(x){
                         fullspCovs <- x[start.cells]
                         if(inherits(x,"RasterLayer")){
                           return(fullspCovs)
                         } else {
                           spCov <- numeric(length(start.cells))
                           zname <- names(attributes(x)$z)
                           if((!zname %in% covNames) & zname!="time") stop("z-value name '",zname,"' must be included in covNames")
                           zvalues <- raster::getZ(x)
                           for(ii in 1:length(zvalues)){
                             spCov[which(rep(moveData[[zname]],each=directions)==zvalues[ii])] <- fullspCovs[which(rep(moveData[[zname]],each=directions)==zvalues[ii]),ii]
                           }
                           return(spCov)
                         }
                       }))
}

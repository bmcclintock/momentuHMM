#' Plot observations on satellite image
#'
#' Plot tracking data on a satellite map. This function plots coordinates in longitude
#' and latitude (not UTM), so if \code{data} coordinates are not provided in longitude and latitude, then the coordinate reference system must be provided using the \code{projargs} argument. This function uses the package \code{ggmap}
#' to fetch a satellite image from Google. An Internet connection is required to use
#' this function.
#'
#' @param data Data frame or \code{\link{momentuHMMData}} object, with necessary fields 'x' (longitudinal direction) and 'y' (latitudinal direction).  A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object is also permitted, from which the data will be extracted.
#' If \code{states=NULL} and a \code{momentuHMM}, \code{miHMM}, or \code{miSum} object is provided, the decoded states are automatically plotted.
#' @param zoom The zoom level, as defined for \code{\link{get_map}}. Integer value between
#' 3 (continent) and 21 (building).
#' @param location Location of the center of the map to be plotted (this must be in the same coordinate reference system as \code{data}).
#' @param segments \code{TRUE} if segments should be plotted between the observations (default),
#' \code{FALSE} otherwise.
#' @param compact \code{FALSE} if tracks should be plotted separately, \code{TRUE}
#' otherwise (default).
#' @param col Palette of colours to use for the dots and segments. If not specified, uses default palette.
#' @param alpha Transparency argument for \code{\link{geom_point}}.
#' @param size Size argument for \code{\link{geom_point}}.
#' @param shape Shape argument for \code{\link{geom_point}}. If \code{states} is provided, then \code{shape} must either be a scalar or a vector of length \code{length(unique(states))}.
#' If \code{states=NULL}, then \code{shape} must either be a scalar or a vector consisting of a value for each individual to be plotted.
#' @param states A sequence of integers, corresponding to the decoded states for these data
#' (such that the observations are colored by states).
#' @param animals Vector of indices or IDs of animals/tracks to be plotted.
#' Default: \code{NULL}; all animals are plotted.
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param return If \code{TRUE}, the function returns a ggplot object (which can be edited and
#' plotted manually). If \code{FALSE}, the function automatically plots the map (default).
#' @param stateNames Optional character vector of length \code{max(states)} indicating state names. Ignored unless \code{states} is provided.
#' @param projargs A character string of PROJ.4 projection arguments indicating the coordinate reference system for \code{data} and \code{location} coordinates (if not longitude and latitude). A \code{\link[sp]{CRS-class}} object is also permitted. If \code{projargs} is provided, the coordinates will be internally transformed to longitude and latitude for plotting.
#'
#' @details If the plot displays the message "Sorry, we have no imagery here", try a
#' lower level of zoom.
#'
#' @references
#' D. Kahle and H. Wickham. ggmap: Spatial Visualization with ggplot2.
#' The R Journal, 5(1), 144-161.
#' URL: http://journal.r-project.org/archive/2013-1/kahle-wickham.pdf
#'
#' @importFrom ggmap get_map ggmap
#' @importFrom ggplot2 geom_point geom_path aes guides scale_color_manual scale_shape_manual
#' @importFrom sp coordinates proj4string spTransform CRS
#' @export

plotSat <- function(data,zoom=NULL,location=NULL,segments=TRUE,compact=TRUE,col=NULL,alpha=1,size=1,shape=16,
                    states=NULL,animals=NULL,ask=TRUE,return=FALSE,stateNames=NULL,projargs=NULL)
{
  #####################
  ## Check arguments ##
  #####################
  
  if(is.null(states) & inherits(data,c("momentuHMM","miHMM","miSum"))){
    if(inherits(data,"momentuHMM")) states <- viterbi(data)
    else if(inherits(data,"miHMM")) states <- data$miSum$Par$states
    else states <- data$Par$states
  }
  if(is.null(stateNames)){
    if(inherits(data,c("momentuHMM","miHMM","miSum"))){
      if(inherits(data,c("momentuHMM","miSum"))) stateNames <- data$stateNames
      else if(inherits(data,"miHMM")) stateNames <- data$miSum$stateNames
    }
  }
  
  if(inherits(data,"momentuHMM")) data <- data$data
  else if(inherits(data,"miHMM")) data <- data$miSum$data
  else if(inherits(data,"miSum")) data <- data$data
  else if(inherits(data,"HMMfits")) stop("data must be a data frame, momentuHMMData, momentuHMM, miHMM, or miSum object")
  
  if(inherits(data,"momentuHMMData")) data <- as.data.frame(data)

  if(is.null(data$x) | is.null(data$y))
    stop("Data should have fields data$x and data$y.")
  
  if(is.null(data$ID))
    data$ID <- rep("Animal1",nrow(data))
  
  if(!is.null(projargs)){
    if(!inherits(projargs,"CRS")) projargs <- sp::CRS(projargs)
    sp::coordinates(data) <- c("x","y")
    sp::proj4string(data) <- projargs
    data <- as.data.frame(sp::spTransform(data,CRS("+proj=longlat +datum=WGS84")))
    if(!is.null(location)){
      if(length(location)!=2)
        stop("'location' should be a vector of two values.")
      location <- data.frame(x=location[1],y=location[2])
      sp::coordinates(location) <- c("x","y")
      sp::proj4string(location) <- projargs
      location <- as.matrix(as.data.frame(sp::spTransform(location,CRS("+proj=longlat +datum=WGS84"))))
    }
  }
  
  if(min(data$x,na.rm=TRUE) < -180 | max(data$x,na.rm=TRUE) > 180 |
     min(data$y,na.rm=TRUE) < -90 | max(data$y,na.rm=TRUE) > 90)
    stop("Coordinates should be longitude/latitude values.")
  
  if(!is.null(zoom)) {
    if(zoom<3 | zoom>21)
      stop("'zoom' should be between 3 (continent) and 21 (building).")
  } else
    stop("'zoom' needs to be specified as an integer between 3 (continent) and 21 (building)")
  
  if(!is.null(location)) {
    if(length(location)!=2)
      stop("'location' should be a vector of two values (longitude,latitude).")
    if(location[1] < -180 | location[1] > 180 | location[2] < -90 | location[2] > 90)
      stop("'location' should be a vector of a longitude value and a latitude value.")
    
    midLon <- location[1]
    midLat <- location[2]
  }
  
  if(length(alpha)!=1)
    stop("'alpha' should be of length 1")
  
  if(!compact & return)
    stop("Cannot return map if not compact. Either set 'compact=TRUE' or 'return=FALSE'.")
  
  x <- y <- NULL # gets rid of no visible binding for global variable 'x' and 'y' NOTE in R cmd check
  ID <- NULL # same for ID
  
  if(is.null(states) & !is.null(stateNames)) stop("'stateNames' cannot be specified unless 'states' is also specified")
  if(!is.null(states)){
    if(length(states)!=nrow(data))
      stop("'states' should have the same length as 'data'")
    if(is.null(stateNames)) stateNames <- sort(unique(states))
    else stateNames <- stateNames[sort(unique(states))]
    data$states <- states
  }
  
  # subset data to tracks to be plotted
  if(!is.null(animals)) {
    if(is.character(animals)) { # animals' IDs provided
      if(!all(animals%in%unique(data$ID))) # ID not found
        stop("Check animals argument.")
      
      data <- subset(data, ID%in%animals)
    }
    if(is.numeric(animals)) { # animals' indices provided
      if(any(animals<1) | any(animals>length(unique(data$ID)))) # index out of bounds
        stop("Check animals argument.")
      
      data <- subset(data, ID%in%unique(ID)[animals])
      animals <- unique(data$ID)
    }
  }
  
  ##############################
  ## Prepare data and colours ##
  ##############################
  # add column for colours
  if(!is.null(states)) { # use "states"
    
    states <- data$states

    nbCol <- length(unique(states))
    data <- cbind(data,col=as.factor(states))
    colname <- "state"
    
    if(length(shape)>1){
      if(length(shape)!=nbCol)
        stop("'shape' should be a scalar or a vector of length ",nbCol)
      data <- cbind(data,shape=as.factor(shape[states]))
    } else {
      data <- cbind(data,shape=as.factor(shape))
      shape <- rep(shape,nbCol)
    }

    
  } else if(compact) { # use "ID"
    if(!is.null(animals))
      nbCol <- length(animals)
    else {
      nbCol <- length(unique(data$ID))
      animals <- unique(data$ID)
    }
    data <- cbind(data,col=as.factor(data$ID))
    colname <- "ID"
    
    if(length(shape)>1){
      if(length(shape)!=nbCol)
        stop("'shape' should be a scalar or a vector of length ",nbCol)
      data <- cbind(data,shape=as.factor(rep(shape,times=table(data$ID))))
    } else {
      data <- cbind(data,shape=as.factor(shape))
      shape <- rep(shape,nbCol)
    }
    
  } else { # no colours needed
    nbCol <- 1
    if(length(shape)>1) stop("'shape' must be a scalar")
    data <- cbind(data,col=as.factor(1),shape=as.factor(shape))
    colname <- ""
  }
  
  # prepare palette
  if(is.null(col)) {
    if(nbCol==1) {
      pal <- "black"
    } else if(nbCol<8) {
      pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    } else {
      # to make sure that all colours are distinct (emulate ggplot default palette)
      hues <- seq(15, 375, length = nbCol + 1)
      pal <- hcl(h = hues, l = 65, c = 100)[1:nbCol]
    }
  } else {
    # if one color given, duplicate for all tracks or states
    if(length(col)==1) {
      col <- rep(col,nbCol)
      nbCol <- 1 # to remove legend below
    }
    
    if(length(col)<nbCol)
      stop(paste("'col' should be at least of length the number of colors needed (",
                 nbCol,").",sep=""))
    
    pal <- col
  }
  
  #################
  ## Plot tracks ##
  #################
  par(ask=ask)
  if(!compact) {
    # loop over tracks
    for(id in unique(data$ID)) {
      subData <- subset(data, ID==id)
      
      # define center of map
      if(is.null(location)) {
        midLon <- (max(subData$x,na.rm=TRUE) + min(subData$x,na.rm=TRUE))/2
        midLat <- (max(subData$y,na.rm=TRUE) + min(subData$y,na.rm=TRUE))/2
      }
      
      # get map of the study region (from google)
      map <- get_map(location=c(lon=midLon, lat=midLat), zoom=zoom, maptype="satellite",
                     source="google")
      
      mapMove <- ggmap(map) + geom_point(aes(x,y,col=col,shape=col),subData,size=size,alpha=alpha)
      
      if(segments)
        mapMove <- mapMove + geom_path(aes(x,y,col=col,group=ID),subData,alpha=alpha)
      
      if(nbCol==1) # no legend if only one colour
        mapMove <- mapMove + scale_color_manual(values=pal) + scale_shape_manual(values=shape) + guides(col=FALSE,shape=FALSE)
      else if(!is.null(stateNames)){
        stateInd <- sort(unique(subData$states))
        mapMove <- mapMove + scale_color_manual(name = colname,labels = stateNames[stateInd],values = pal[stateInd]) + scale_shape_manual(name = colname,labels = stateNames[stateInd],values = shape[stateInd])
      } else {
        mapMove <- mapMove + scale_color_manual(name = colname,labels = animals,values = pal) + scale_shape_manual(name = colname,labels = animals,values = shape)
      }
      plot(mapMove) # plot map
    }
  } else {
    # define center of map
    if(is.null(location)) {
      midLon <- (max(data$x,na.rm=TRUE) + min(data$x,na.rm=TRUE))/2
      midLat <- (max(data$y,na.rm=TRUE) + min(data$y,na.rm=TRUE))/2
    }
    
    # get map of the study region (from google)
    map <- get_map(location=c(lon=midLon, lat=midLat), zoom=zoom, maptype="satellite",
                   source="google")
    
    mapMove <- ggmap(map) + geom_point(aes(x,y,col=col,shape=col),data,size=size,alpha=alpha)
    
    if(segments)
      mapMove <- mapMove + geom_path(aes(x,y,col=col,group=ID),data,alpha=alpha)
    
    if(nbCol==1) # no legend if only one colour
      mapMove <- mapMove + scale_color_manual(values=pal) + scale_shape_manual(values=shape) + guides(col=FALSE,shape=FALSE)
      else if(!is.null(stateNames)){
        mapMove <- mapMove + scale_color_manual(name = colname,labels = stateNames,values = pal) + scale_shape_manual(name = colname,labels = stateNames,values = shape)
      } else {
        mapMove <- mapMove + scale_color_manual(name = colname,labels = animals,values = pal) + scale_shape_manual(name = colname,labels = animals,values = shape)
      }
    if(return) {
      par(ask=FALSE)
      return(mapMove) # return map object
    } else
      plot(mapMove) # plot map
  }
  
  par(ask=FALSE)
}

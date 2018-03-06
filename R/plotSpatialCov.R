#' Plot observations on raster image
#'
#' Plot tracking data over a raster layer.
#'
#' @param data Data frame of the location data, with necessary fields 'x' (longitudinal direction) and
#' 'y' (latitudinal direction).
#' @param spatialCov \code{\link[raster]{RasterLayer-class}} object on which to plot the location data
#' @param segments \code{TRUE} if segments should be plotted between the observations (default),
#' \code{FALSE} otherwise.
#' @param compact \code{FALSE} if tracks should be plotted separately, \code{TRUE}
#' otherwise (default).
#' @param col Palette of colours to use for the dots and segments. If not specified, uses default palette.
#' @param alpha Transparency argument for \code{\link{geom_point}}.
#' @param size Size argument for \code{\link{geom_point}}.
#' @param states A sequence of integers, corresponding to the decoded states for these data. If
#' specified, the observations are colored by states.
#' @param animals Vector of indices or IDs of animals/tracks to be plotted.
#' Default: \code{NULL}; all animals are plotted.
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param return If \code{TRUE}, the function returns a ggplot object (which can be edited and
#' plotted manually). If \code{FALSE}, the function automatically plots the map (default).
#' 
#' @examples 
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' 
#' # plot simulated data over forest raster automatically loaded with the packge
#' spatialCov<-list(forest=forest)
#' data <- simData(nbAnimals=2,nbStates=2,dist=list(step=stepDist,angle=angleDist),
#'                 Par=list(step=c(100,1000,50,100),angle=c(0,0,0.1,5)),
#'                 beta=matrix(c(5,-10,-25,50),nrow=2,ncol=2,byrow=TRUE),
#'                 formula=~forest,spatialCovs=spatialCov,
#'                 obsPerAnimal=225,states=TRUE)
#' 
#' plotSpatialCov(data,forest,states=data$states)
#'
#' @importFrom ggplot2 ggplot geom_point geom_path aes geom_raster guides guide_legend theme
#' @importFrom ggplot2 element_rect element_blank scale_color_manual coord_equal labs
#' @importFrom raster rasterToPoints
#' @export

plotSpatialCov <- function(data,spatialCov,segments=TRUE,compact=TRUE,col=NULL,alpha=1,size=1,
                           states=NULL,animals=NULL,ask=TRUE,return=FALSE)
{
  #####################
  ## Check arguments ##
  #####################
  if(is.null(data$x) | is.null(data$y))
    stop("Data should have fields data$x and data$y.")
  
  if(length(alpha)!=1)
    stop("'alpha' should be of length 1")
  
  if(!compact & return)
    stop("Cannot return map if not compact. Either set 'compact=TRUE' or 'return=FALSE'.")
  
  if(!inherits(spatialCov,"RasterLayer")) 
    stop("spatialCov should be of class 'RasterLayer'")
  
  x <- y <- NULL #gets rid of no visible binding for global variable 'x' and 'y' NOTE in R cmd check
  ID <- NULL # same for ID
  
  #############################
  ## Prepare data and colors ##
  #############################
  # add column for colours
  if(!is.null(states)) {
    # color by "states"
    
    if(length(states)!=nrow(data))
      stop("'states' should have the same length as 'data'")
    
    nbCol <- length(unique(states))
    data <- cbind(data,col=as.factor(states))
    colname <- "state"
    
  } else if(compact & !is.null(data$ID)) {
    # color by "ID"
    
    if(!is.null(animals))
      nbCol <- length(animals)
    else
      nbCol <- length(unique(data$ID))
    data <- cbind(data,col=as.factor(data$ID))
    colname <- "ID"
    
  } else {
    # no colors needed
    
    if(is.null(data$ID)) data$ID <- rep(1,nrow(data))
    
    nbCol <- 1
    data <- cbind(data,col=as.factor(data$ID))
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
    }
  }
  
  #################
  ## Plot tracks ##
  #################
  par(ask=ask)
  
  map.cov <- raster::rasterToPoints(spatialCov)
  dfcov <- data.frame(map.cov)
  colnames(dfcov) <- c("x", "y", names(spatialCov))
  
  if(!compact) {
    #loop over tracks
    for(id in unique(data$ID)) {
      subData <- subset(data,ID==id)
      
      mapMove <- ggplot(dfcov, aes(x=x, y=y)) + 
        theme(panel.background = element_rect(fill=NA)) +
        theme(panel.border = element_rect(colour = "black",fill=NA)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_raster(data = dfcov,aes(fill=dfcov[[names(spatialCov)]])) +
        labs(fill = names(spatialCov)) +
        geom_point(aes(x,y,col=col),subData,alpha=alpha,size=size) +
        coord_equal()
      
      if(segments)
        mapMove <- mapMove + geom_path(aes(x,y,col=col,group=ID),subData,alpha=alpha)
      
      if(nbCol==1) # no legend if only one color
        mapMove <- mapMove + scale_color_manual(values=pal) + guides(col=FALSE)
      else
        mapMove <- mapMove + scale_color_manual(values=pal, name=colname)
      
      plot(mapMove)
    }
  } else {
    mapMove <- ggplot(dfcov, aes(x=x, y=y)) + 
      theme(panel.background = element_rect(fill=NA)) +
      theme(panel.border = element_rect(colour = "black",fill=NA)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_raster(data = dfcov,aes(fill=dfcov[[names(spatialCov)]])) +
      labs(fill = names(spatialCov)) +
      geom_point(aes(x,y,col=col),data,alpha=alpha,size=size) +
      coord_equal()
    
    if(segments)
      mapMove <- mapMove + geom_path(aes(x,y,col=col,group=ID),data,alpha=alpha)
    
    if(nbCol==1) # no legend if only one color
      mapMove <- mapMove + scale_color_manual(values=pal) + guides(col=FALSE)
    else
      mapMove <- mapMove + scale_color_manual(values=pal, name=colname)
    
    if(return) {
      par(ask=FALSE)
      return(mapMove) # return map object
    } else
      plot(mapMove) # plot map
  }
  
  par(ask=FALSE)
}


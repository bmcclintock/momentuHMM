#' Plot observations on raster image
#'
#' Plot tracking data over a raster layer.
#'
#' @param data Data frame or \code{\link{momentuHMMData}} object, with necessary fields 'x' (longitudinal direction) and 'y' (latitudinal direction).  A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object is also permitted, from which the data will be extracted.
#' If \code{states=NULL} and a \code{momentuHMM}, \code{miHMM}, or \code{miSum} object is provided, the decoded states are automatically plotted.
#' @param spatialCov \code{\link[=RasterLayer-class]{RasterLayer}} object on which to plot the location data
#' @param segments \code{TRUE} if segments should be plotted between the observations (default),
#' \code{FALSE} otherwise.
#' @param compact \code{FALSE} if tracks should be plotted separately, \code{TRUE}
#' otherwise (default).
#' @param col Palette of colours to use for the dots and segments. If not specified, uses default palette.
#' @param alpha Transparency argument for \code{\link{geom_point}}.
#' @param size Size argument for \code{\link{geom_point}}.
#' @param shape Shape argument for \code{\link{geom_point}}. If \code{states} is provided, then \code{shape} must either be a scalar or a vector of length \code{length(unique(states))}.
#' If \code{states=NULL}, then \code{shape} must either be a scalar or a vector consisting of a value for each individual to be plotted.
#' @param states A sequence of integers, corresponding to the decoded states for these data. If
#' specified, the observations are colored by states.
#' @param animals Vector of indices or IDs of animals/tracks to be plotted.
#' Default: \code{NULL}; all animals are plotted.
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param return If \code{TRUE}, the function returns a ggplot object (which can be edited and
#' plotted manually). If \code{FALSE}, the function automatically plots the map (default).
#' @param stateNames Optional character vector of length \code{max(states)} indicating state names. Ignored unless \code{states} is provided.
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
#' @importFrom ggplot2 element_rect element_blank scale_color_manual scale_shape_manual coord_equal labs
#' @importFrom raster rasterToPoints
#' @export

plotSpatialCov <- function(data,spatialCov,segments=TRUE,compact=TRUE,col=NULL,alpha=1,size=1,shape=16,
                           states=NULL,animals=NULL,ask=TRUE,return=FALSE,stateNames=NULL)
{
  #####################
  ## Check arguments ##
  #####################
  if(!inherits(spatialCov,"RasterLayer")) 
    stop("spatialCov should be of class 'RasterLayer'")
  
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
  
  if(inherits(data,"momentuHMMData") | inherits(data,"momentuHierHMMData")) data <- as.data.frame(data)
  
  coordNames <- c("x","y")
  if(!is.null(attr(data,'coords'))) {
    coordNames <- attr(data,'coords')
  }
  if(is.null(data[[coordNames[1]]]) | is.null(data[[coordNames[2]]]))
    stop("Data should have fields data$",coordNames[1]," and data$",coordNames[2],".")
  
  data$x <- data[[coordNames[1]]]
  data$y <- data[[coordNames[2]]]
  
  if(length(alpha)!=1)
    stop("'alpha' should be of length 1")
  
  if(!compact & return)
    stop("Cannot return plot if not compact. Either set 'compact=TRUE' or 'return=FALSE'.")
  
  if(is.null(data$ID))
    data$ID <- rep("Animal1",nrow(data)) # default ID if none provided
  
  x <- y <- NULL #gets rid of no visible binding for global variable 'x' and 'y' NOTE in R cmd check
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
      pal <- "#E69F00"
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
  map.cov <- raster::rasterToPoints(spatialCov)
  dfcov <- data.frame(map.cov)
  colnames(dfcov) <- c("x", "y", names(spatialCov))
  
  if(!compact) {
    
    par(ask=ask)
    
    #loop over tracks
    for(id in unique(data$ID)) {
      subData <- subset(data,ID==id)
      
      mapMove <- ggplot(dfcov, aes(x=x, y=y)) + 
        theme(panel.background = element_rect(fill=NA)) +
        theme(panel.border = element_rect(colour = "black",fill=NA)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_raster(data = dfcov,aes(fill=dfcov[[names(spatialCov)]])) +
        labs(fill = names(spatialCov)) +
        geom_point(aes(x,y,col=col,shape=col),subData,alpha=alpha,size=size,show.legend=ifelse(nbCol>1,TRUE,FALSE)) +
        coord_equal()
      
      if(segments)
        mapMove <- mapMove + geom_path(aes(x,y,col=col,group=ID),subData,alpha=alpha)
      
      if(nbCol==1) # no legend if only one color
        mapMove <- mapMove + scale_color_manual(values=pal) + scale_shape_manual(values=shape) + guides(col=FALSE)
      else if(!is.null(stateNames)){
        stateInd <- sort(unique(subData$states))
        mapMove <- mapMove + scale_color_manual(name = colname,labels = stateNames[stateInd],values = pal[stateInd]) + scale_shape_manual(name = colname,labels = stateNames[stateInd],values = shape[stateInd])
      } else {
        mapMove <- mapMove + scale_color_manual(name = colname,labels = animals,values = pal) + scale_shape_manual(name = colname,labels = animals,values = shape)
      }
      
      plot(mapMove)
    }
  } else {
    mapMove <- ggplot(dfcov, aes(x=x, y=y)) + 
      theme(panel.background = element_rect(fill=NA)) +
      theme(panel.border = element_rect(colour = "black",fill=NA)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_raster(data = dfcov,aes(fill=dfcov[[names(spatialCov)]])) +
      labs(fill = names(spatialCov)) +
      geom_point(aes(x,y,col=col,shape=col),data,alpha=alpha,size=size,show.legend=ifelse(nbCol>1,TRUE,FALSE)) +
      coord_equal()
    
    if(segments)
      mapMove <- mapMove + geom_path(aes(x,y,col=col,group=ID),data,alpha=alpha)
    
    if(nbCol==1) # no legend if only one color
      mapMove <- mapMove + scale_color_manual(values=pal) + scale_shape_manual(values=shape) + guides(col=FALSE)
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


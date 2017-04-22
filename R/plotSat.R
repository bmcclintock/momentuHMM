#' Plot observations on satellite image
#'
#' Plot tracking data on a satellite map. This function only works with longitude
#' and latitude values (not with UTM coordinates), and uses the package \code{ggmap}
#' to fetch a satellite image from Google. An Internet connection is required to use
#' this function.
#'
#' @param data Data frame of the data, with necessary fields 'x' (longitude values) and
#' 'y' (latitude values).
#' @param zoom The zoom level, as defined for \code{\link{get_map}}. Integer value between
#' 3 (continent) and 21 (building).
#' @param location Location of the center of the map to be plotted.
#' @param segments \code{TRUE} if segments should be plotted between the observations (default),
#' \code{FALSE} otherwise.
#' @param compact \code{FALSE} if tracks should be plotted separately, \code{TRUE}
#' otherwise (default).
#' @param col Palette of colours to use for the dots and segments. If not specified, uses default palette.
#' @param alpha Transparency argument for \code{\link{geom_point}}.
#' @param size Size argument for \code{\link{geom_point}}.
#' @param states A sequence of integers, corresponding to the decoded states for these data
#' (such that the observations are colored by states).
#' @param animals Vector of indices or IDs of animals/tracks to be plotted.
#' Default: \code{NULL}; all animals are plotted.
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param return If \code{TRUE}, the function returns a ggplot object (which can be edited and
#' plotted manually). If \code{FALSE}, the function automatically plots the map (default).
#' @param stateNames Optional character vector of length \code{unique(states)} indicating state names. Ignored unless \code{states} is provided.
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
#' @importFrom ggplot2 geom_point geom_path aes guides scale_color_manual
#' @export

plotSat <- function(data,zoom=NULL,location=NULL,segments=TRUE,compact=TRUE,col=NULL,alpha=1,size=1,
                    states=NULL,animals=NULL,ask=TRUE,return=FALSE,stateNames=NULL)
{
  #####################
  ## Check arguments ##
  #####################
  if(is.null(data$x) | is.null(data$y))
    stop("Data should have fields data$x and data$y.")
  
  if(is.null(data$ID))
    data$ID <- rep("Animal1",nrow(data))
  
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
  
  ##############################
  ## Prepare data and colours ##
  ##############################
  # add column for colours
  if(!is.null(states)) { # use "states"
    if(length(states)!=nrow(data))
      stop("'states' should have the same length as 'data'")
    
    nbCol <- length(unique(states))
    data <- cbind(data,col=as.factor(states))
    colname <- "state"
    
  } else if(compact) { # use "ID"
    if(!is.null(animals))
      nbCol <- length(animals)
    else
      nbCol <- length(unique(data$ID))
    data <- cbind(data,col=as.factor(data$ID))
    colname <- "ID"
    
  } else { # no colours needed
    nbCol <- 1
    data <- cbind(data,col=as.factor(1))
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
      
      mapMove <- ggmap(map) + geom_point(aes(x,y,col=col),subData,size=size,alpha=alpha)
      
      if(segments)
        mapMove <- mapMove + geom_path(aes(x,y,col=col),subData,alpha=alpha)
      
      if(nbCol==1) # no legend if only one colour
        mapMove <- mapMove + scale_color_manual(values=pal) + guides(col=FALSE)
      else
        mapMove <- mapMove + scale_color_manual(name = colname, values=pal)
      
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
    
    mapMove <- ggmap(map) + geom_point(aes(x,y,col=col),data,size=size,alpha=alpha)
    
    if(segments)
      mapMove <- mapMove + geom_path(aes(x,y,col=col,group=ID),data,alpha=alpha)
    
    if(nbCol==1) # no legend if only one colour
      mapMove <- mapMove + scale_color_manual(values=pal) + guides(col=FALSE)
    else {
      if(!is.null(stateNames)) mapMove <- mapMove + scale_color_manual(name = colname, values=pal, labels=stateNames)
      else mapMove <- mapMove + scale_color_manual(name = colname, values=pal)
    }
    if(return) {
      par(ask=FALSE)
      return(mapMove) # return map object
    } else
      plot(mapMove) # plot map
  }
  
  par(ask=FALSE)
}

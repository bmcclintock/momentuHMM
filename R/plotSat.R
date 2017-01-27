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
#' @param col Color(s) of the dots and segments. Should be either of length 1, or of
#' the length of the data.
#' @param alpha Transparency argument for \code{\link{geom_point}}.
#' @param size Size argument for \code{\link{geom_point}}.
#' @param states A sequence of integers, corresponding to the decoded states for these data
#' (such that the observations are colored by states). If 'states' if specified, the
#' argument 'col' gives the colors correponding to each state.
#' @param animals Vector of indices or IDs of animals/tracks to be plotted.
#' Default: \code{NULL}; all animals are plotted.
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param return If \code{TRUE}, the function returns a ggplot object (which can be edited and
#' plotted manually). If \code{FALSE}, the function automatically plots the map (default).
#'
#' @details If the plot displays the message "Sorry, we have no imagery here", try a
#' lower level of zoom.
#'
#' @references
#' D. Kahle and H. Wickham. ggmap: Spatial Visualization with ggplot2.
#' The R Journal, 5(1), 144-161.
#' URL: http://journal.r-project.org/archive/2013-1/kahle-wickham.pdf
#'
#' @importFrom ggmap get_map
#' @importFrom ggmap ggmap
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 aes
#' @export

plotSat <- function(data,zoom=NULL,location=NULL,segments=TRUE,compact=TRUE,col=1,alpha=1,size=1,
                    states=NULL,animals=NULL,ask=TRUE,return=FALSE)
{
  #####################
  ## Check arguments ##
  #####################
  if(is.null(data$x) | is.null(data$y))
    stop("Data should have fields data$x and data$y.")
  
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
  }
  
  if(is.null(states)) {
    if(length(col)!=1 & length(col)!=nrow(data))
      stop("'col' should be of length 1, or length of the data")
  } else {
    if(length(states)!=nrow(data))
      stop("'states' should have the same length as 'data'")
    
    # number of states
    nbStates <- length(unique(states))
    
    if(length(col)!=nbStates) {
      if(nbStates<8) {
        pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        col <- pal[states]
      } else
        col <- rainbow(nbStates)[states] # to make sure that all colors are distinct
    } else {
      col <- col[states]
    }
  }
  
  if(length(alpha)!=1)
    stop("'alpha' should be of length 1")
  
  #################
  ## Plot tracks ##
  #################
  # color by track if several tracks plotted on a single map
  # (does not use col==1 to avoid "condition has length > 1" warning message)
  if(is.null(states) & compact & !is.null(data$ID) & length(col)==1 & col[1]==1) {
    if(length(unique(data$ID))<8)
      pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    else
      pal <- rainbow(length(unique(data$ID)))
    
    col <- rep(NA,nrow(data))
    for(i in 1:length(unique(data$ID)))
      col[which(data$ID==unique(data$ID)[i])] <- pal[i]
  }
  
  # is there one color value per point? used later
  manyCols <- length(col)==nrow(data)
  
  # select subset of data to be plotted
  if(!is.null(animals)) {
    if(is.character(animals)) { # animals' IDs provided
      if(manyCols)
        col <- col[which(data$ID%in%animals)]
      
      data <- data[which(data$ID%in%animals),]
    }
    
    
    if(is.numeric(animals)) { # animals' indices provided
      nbAnimals <- length(unique(data$ID))
      if(length(which(animals<1))>0 | length(which(animals>nbAnimals))>0) # index out of bounds
        stop("Check 'animals' argument, index out of bounds")
      
      if(manyCols)
        col <- col[which(data$ID%in%unique(data$ID)[animals])]
      
      data <- data[which(data$ID%in%unique(data$ID)[animals]),]
    }
  }
  
  # define number of tracks to be plotted
  if(is.null(data$ID) | compact)
    nbTracks <- 1
  else
    nbTracks <- length(unique(data$ID))
  
  # loop over tracks
  for(zoo in 1:nbTracks) {
    if(nbTracks==1) {
      subData <- data
      subCol <- col
    } else {
      par(ask=ask)
      
      # select data for track 'zoo'
      subData <- data[which(data$ID==unique(data$ID)[zoo]),]
      
      # select colors for track 'zoo'
      if(manyCols)
        subCol <- col[which(data$ID==unique(data$ID)[zoo])]
      else
        subCol <- col
    }
    
    # define center of map
    if(is.null(location)) {
      midLon <- min(subData$x,na.rm=TRUE) + (max(subData$x,na.rm=TRUE)-min(subData$x,na.rm=TRUE))/2
      midLat <- min(subData$y,na.rm=TRUE) + (max(subData$y,na.rm=TRUE)-min(subData$y,na.rm=TRUE))/2
    } else {
      midLon <- location[1]
      midLat <- location[2]
    }
    
    # get map of the study region (from google)
    map <- get_map(location=c(lon=midLon, lat=midLat), zoom=zoom, maptype="satellite",
                   source="google")
    
    # define map with dots
    mapMove <- ggmap(map) + geom_point(aes_string(x="x", y="y"), data=subData, col=subCol,
                                       alpha=alpha,size=size)
    
    if(segments) {
      # if several tracks on one map
      if(compact & !is.null(data$ID)) {
        # loop over tracks
        for(id in unique(data$ID)) {
          ind <- which(data$ID==id)
          
          xto <- subData$x[ind[-1]]
          yto <- subData$y[ind[-1]]
          xfrom <- subData$x[ind[-length(ind)]]
          yfrom <- subData$y[ind[-length(ind)]]
          
          if(manyCols)
            subCol <- col[ind[-length(ind)]]
          
          seg <- data.frame(xfrom=xfrom,yfrom=yfrom,xto=xto,yto=yto)
          
          # add segments to map
          mapMove <- mapMove + geom_segment(aes(x=xfrom, y=yfrom, xend=xto, yend=yto), data=seg,
                                            col=subCol, alpha=alpha)
        }
      } else {
        # define segments to plot
        xto <- subData$x[-1]
        yto <- subData$y[-1]
        xfrom <- subData$x[-nrow(subData)]
        yfrom <- subData$y[-nrow(subData)]
        
        # remove last color, because there is one fewer segment
        if(manyCols)
          subCol <- subCol[-length(subCol)]
        
        seg <- data.frame(xfrom=xfrom,yfrom=yfrom,xto=xto,yto=yto)
        
        # add segments to map
        mapMove <- mapMove + geom_segment(aes(x=xfrom, y=yfrom, xend=xto, yend=yto), data=seg,
                                          col=subCol, alpha=alpha)
      }
    }
    
    if(return)
      return(mapMove) # return map object
    else
      plot(mapMove) # plot map
  }
  
  par(ask=FALSE)
}
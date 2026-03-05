
#' Convert standard direction angles (in radians relative to the x-axis) to turning angle covariates suitable for circular-circular regression on the angle mean
#' 
#' This function can be used to convert angular covariates (e.g., ocean currents, wind direction) measured in radians relative to the x-axis to turning angle
#' covariates sutiable for circular-circular regression in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}.
#' 
#' @param refAngle Numeric vector of standard direction angles (in radians) relative to the x-axis, where 0 = east, pi/2 = north, pi = west, -pi/2 = south
#' @param data data frame containing fields for the x- and y-coordinates (identified by \code{coordNames}) and 'ID' (if more than one individual)
#' @param coordNames Names of the columns of coordinates in \code{data}. Default: \code{c("x","y")}.
#' 
#' @return A vector of turning angles between the movement direction at time step t-1 and \code{refAngle} at time t
#' 
#' @examples
#' # extract data from momentuHMM example
#' data<-example$m$data
#' 
#' # generate fake angle covariates
#' u <- rnorm(nrow(data)) # horizontal component
#' v <- rnorm(nrow(data)) # vertical component
#' refAngle <- atan2(v,u)
#' 
#' # add turning angle covariate to data
#' data$cov3 <- circAngles(refAngle=refAngle,data=data)
#' 
#' @export

circAngles <- function(refAngle, data, coordNames = c("x", "y")) {
  
  x <- data[[coordNames[1]]]
  y <- data[[coordNames[2]]]
  n <- length(x)
  
  ID <- data$ID
  if (is.null(ID)) ID <- rep("Animal1", n)
  id_firsts <- match(unique(ID), ID)
  
  x_prev <- c(x[1], x[-n])
  y_prev <- c(y[1], y[-n])
  
  x_prev[id_firsts] <- x[id_firsts]
  y_prev[id_firsts] <- y[id_firsts]
  
  #if (type == "LL") {
  #  rad <- pi / 180
  #  lon1 <- x_prev * rad
  #  lat1 <- y_prev * rad
  #  lon2 <- x * rad
  #  lat2 <- y * rad
  #  
  #  dLon <- lon2 - lon1
  #  by <- sin(dLon) * cos(lat2)
  #  bx <- cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon)
  #  heading <- atan2(by, bx)
  #  
  #} else if (type == "UTM") {
    dx <- x - x_prev
    dy <- y - y_prev
    heading <- atan2(dy, dx)
  #} else {
  #  stop("type must be 'LL' or 'UTM'")
  #}
  
  angle <- ((heading - refAngle + pi) %% (2 * pi)) - pi
  
  return(angle)
}

#refAngles<-function(x,y){
#  angle<-atan2(y[1]-y[2],x[1]-x[2])
#  while(angle<=(-pi)) angle <- angle + 2*pi
#  while(angle>pi) angle <- angle -2*pi
#  return(angle)
#}

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

circAngles<-function(refAngle,data,coordNames=c("x","y")){
  
  if(!is.data.frame(data)) stop("data must be a data frame")
  if(any(dim(data)==0)) stop("data is empty")
  if(length(coordNames)!=2) stop('coordNames must be of length 2')
  
  if(!is.null(data$ID)) ID <- as.character(data$ID) # homogenization of numeric and string IDs
  else ID <- rep(1,nrow(data)) # default ID if none provided
  
  if(nrow(data)!=length(refAngle)) stop('refAngle must be of length ',nrow(data))
  
  if(min(refAngle)<= -pi | max(refAngle)> pi) stop('refAngle must be in (-pi,pi]')
  
  x <- data[[coordNames[1]]]
  y <- data[[coordNames[2]]]
  
  ind<-as.numeric(table(ID))
  cumind<-c(0,cumsum(ind))
  angle <- rep(0,length(x))
  for(s in 1:length(ind)){
    for(i in cumind[s]+2:ind[s]){
      w <- c(x[i]-x[i-1],y[i]-y[i-1])
      #b <- -(atan2(w[2],w[1])-pi/2) # bearing; 0 is north, pi/2 east
      #b <- atan2(sin(b),cos(b)) # bearing on [-pi,pi)
      b <- atan2(w[2],w[1]) # 0 is east, pi/2 north, -pi/2 south, pi west
      angle[i] <- refAngle[i] - b
      while(angle[i]<=(-pi)) angle[i] <- angle[i] + 2*pi
      while(angle[i]>pi) angle[i] <- angle[i] -2*pi
    }
  }
  class(angle) <- c(class(angle), "angle")
  return(angle)
}

#refAngles<-function(x,y){
#  angle<-atan2(y[1]-y[2],x[1]-x[2])
#  while(angle<=(-pi)) angle <- angle + 2*pi
#  while(angle>pi) angle <- angle -2*pi
#  return(angle)
#}

#' Calculate distance between points y and z and turning angle between points x, y, and z
#' 
#' @param x location 1
#' @param y location 2
#' @param z location 3
#' @param type \code{'UTM'} if easting/northing provided (the default), \code{'LL'} if longitude/latitude
#' 
#' @return 2-vector with first element the distance between y and z and second element the turning angle
#' between (x,y) and (y,z).
#' 
#' @details 
#' Used in \code{\link{prepData}} and \code{\link{simData}} to get distance and turning angle covariates 
#' between locations (x1,x2), (y1,y2) and activity center (z1,z2).
#' 
#' If \code{type='LL'} then distance is calculated as great circle distance using \code{\link[sp]{spDistsN1}}, and turning angle is calculated based on initial bearings using \code{\link[geosphere]{bearing}}.
#' @importFrom sp spDistsN1

distAngle<-function(x,y,z,type='UTM'){
  dist <- spDistsN1(pts = matrix(y,ncol=2),
                    pt = z,
                    longlat = (type=='LL'))
  angle <- turnAngle(x,y,z,type)
  c(dist,angle)
}
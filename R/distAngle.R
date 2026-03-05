
#' Calculate distance between points y and z and turning angle between points x, y, and z
#' 
#' @param x location 1
#' @param y location 2
#' @param z location 3
#' @param type \code{'UTM'} if easting/northing provided (the default), \code{'LL'} if longitude/latitude
#' @param angleCov logical indicating to not return NA when x=y or y=z. Default: TRUE (i.e. NA is not returned if x=y or y=z).
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

distAngle <- function(x, y, z, type = "UTM", angleCov = TRUE) {
  
  if (is.null(dim(x))) x <- matrix(x, ncol = 2)
  if (is.null(dim(y))) y <- matrix(y, ncol = 2)
  if (is.null(dim(z))) z <- matrix(z, ncol = 2)
  
  n <- nrow(y)
  
  if (nrow(z) == 1 && n > 1) z <- matrix(rep(z, each = n), ncol = 2)
  if (nrow(x) == 1 && n > 1) x <- matrix(rep(x, each = n), ncol = 2)
  
  if (type == "LL") {
    rad <- pi / 180
    
    dist <- sp::spDists(y, z, longlat = TRUE, diagonal = TRUE)
    
    # Great Circle Heading from x to y
    dLon1 <- (y[, 1] - x[, 1]) * rad
    by1 <- sin(dLon1) * cos(y[, 2] * rad)
    bx1 <- cos(x[, 2] * rad) * sin(y[, 2] * rad) - 
      sin(x[, 2] * rad) * cos(y[, 2] * rad) * cos(dLon1)
    phi1 <- atan2(by1, bx1)
    
    # Great Circle Heading from y to z
    dLon2 <- (z[, 1] - y[, 1]) * rad
    by2 <- sin(dLon2) * cos(z[, 2] * rad)
    bx2 <- cos(y[, 2] * rad) * sin(z[, 2] * rad) - 
      sin(y[, 2] * rad) * cos(z[, 2] * rad) * cos(dLon2)
    phi2 <- atan2(by2, bx2)
    
  } else if (type == "UTM") {
    dist <- sqrt((z[, 1] - y[, 1])^2 + (z[, 2] - y[, 2])^2)
    
    phi1 <- atan2(y[, 2] - x[, 2], y[, 1] - x[, 1])
    phi2 <- atan2(z[, 2] - y[, 2], z[, 1] - y[, 1])
  } else {
    stop("type must be 'LL' or 'UTM'")
  }
  
  zero_step <- (x[, 1] == y[, 1] & x[, 2] == y[, 2])
  phi1[zero_step] <- 0
  
  angle <- ((phi2 - phi1 + pi) %% (2 * pi)) - pi
  
  if (angleCov) {
    res <- cbind(dist = dist, angle = angle)
    if (n == 1) return(as.numeric(res)) 
    return(res)
  } else {
    if (n == 1) return(as.numeric(dist))
    return(dist)
  }
}
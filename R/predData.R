#' Create \code{momentuHMMData} object from \code{crwData} object 
#'
#' Create a \code{\link{momentuHMMData}} object from a \code{\link{crwData}} object based on the best predicted 
#' locations (\code{crwData$crwPredict$mu.x} and \code{crwData$crwPredict$mu.y}) and any additional covariates
#' 
#' @param crwData A \code{\link{crwData}} object (as returned by \code{\link{crawlWrap}})
#' @param covNames Names of any covariates in \code{crwData$crwPredict}. See \code{\link{prepData}}. 
#' @param spatialCovs List of raster layer(s) for any spatial covariates not included in \code{crwData$crwPredict}. See \code{\link{prepData}}. 
#' @param centers 2-column matrix providing the x-coordinates (column 1) and y-coordinates (column 2) for any activity centers (e.g., potential 
#' centers of attraction or repulsion) from which distance and angle covariates will be calculated. See \code{\link{prepData}}. 
#' @param angleCovs Character vector indicating the names of any circular-circular regression angular covariates in \code{crwData$crwPredict} that need conversion from standard direction (in radians relative to the x-axis) to turning angle (relative to previous movement direction) 
#' See \code{\link{prepData}}.
#' 
#' @return A \code{\link{momentuHMMData}} object suitable for \code{\link{fitHMM}}.
#' @export
predData<-function(crwData,covNames=NULL,spatialCovs=NULL,centers=NULL,angleCovs=NULL){
  if(is.crwData(crwData)){
    track <- MIfitHMM(crwData,nSims=1,ncores=1,covNames=covNames,spatialCovs=spatialCovs,centers=centers,angleCovs=angleCovs,fit=FALSE)[[1]]
  } else stop("crwData must be a crwData object")
}
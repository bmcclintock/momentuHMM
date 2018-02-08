
#' Merge crwData object with additional data streams and/or covariates
#' 
#' This function can be used to merge \code{\link{crwData}} objects (as returned by \code{\link{crawlWrap}}) with additional data streams
#' and/or covariates that are unrelated to location. 
#' 
#' Specifically, the function merges the \code{crwData$crwPredict} data frame with \code{data}
#' based on the \code{ID} and \code{Time.name} columns.  Thus both \code{crwData$crwPredict} and \code{data} must contain \code{ID} and \code{Time.name} columns. 
#' 
#' Only rows of \code{data} with \code{ID} and \code{Time.name} values that exactly match \code{crwData$crwPredict} are merged. Typically, the \code{Time.name} column
#' in \code{data} should match predicted times of locations in \code{crwData$crwPredict} (i.e. those corresponding to \code{crwData$crwPredict$locType=="p"})
#' 
#' @param crwData A \code{\link{crwData}} object
#' @param data A data frame containing required columns \code{ID} and \code{Time.name}, plus any additional data streams and/or covariates to merge with \code{crwData}.
#' @param Time.name Character string indicating name of the time column to be used for merging
#' 
#' @return A \code{\link{crwData}} object
#' @examples
#' \dontrun{
#' # extract simulated obsData from example data
#' obsData <- miExample$obsData
#' 
#' # extract crwMLE inputs from example data
#' inits <- miExample$inits # initial state
#' err.model <- miExample$err.model # error ellipse model
#'
#' # Fit crwMLE models to obsData and predict locations 
#' # at default intervals for both individuals
#' crwOut <- crawlWrap(obsData=obsData,
#'          theta=c(4,0),fixPar=c(1,1,NA,NA),
#'          initial.state=inits,
#'          err.model=err.model,attempts=100)
#'          
#' # create data frame with fake data stream
#' data <- data.frame(ID=rep(factor(c(1,2)),times=c(753,652)),
#'                    time=c(1:753,1:652),
#'                    fake=rpois(753+652,5))
#' 
#' # merge fake data stream with crwOut
#' crwOut <- crawlMerge(crwOut,data,"time")
#' }
#' @export
crawlMerge<-function(crwData, data, Time.name){
  if(!is.crwData(crwData)) stop("crwData must be a crwData object")
  crwFits <- crwData$crwFits
  if(!(Time.name %in% names(crwData$crwPredict))) stop(Time.name," not found in crwData$crwPredict")
  if(!("ID" %in% names(data))) stop("ID not found in data")
  if(!(Time.name %in% names(data))) stop(Time.name," not found in data")
  covNames <- names(data)[!(names(data) %in% c("ID",Time.name))]
  if(any(covNames %in% names(crwData$crwPredict)[!(names(crwData$crwPredict) %in% c("ID",Time.name))])) stop("Data stream and/or covariate names in data cannot match column names in crwData$crwPredict")
  tmp<-merge(crwData$crwPredict,data,by=c("ID",Time.name),all.x=TRUE,all.y=FALSE)
  crwData$crwPredict[covNames]<-tmp[covNames]
  return(crwData)
}
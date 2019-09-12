
#' Merge crwData or crwHierData object with additional data streams and/or covariates
#' 
#' This function can be used to merge \code{\link{crwData}} or \code{\link{crwHierData}} objects (as returned by \code{\link{crawlWrap}}) with additional data streams
#' and/or covariates that are unrelated to location. 
#' 
#' Specifically, the function merges the \code{crwData$crwPredict} data frame with \code{data}
#' based on the \code{ID}, \code{Time.name}, and, if \code{crwData} is hierarchical, \code{level} columns.  Thus \code{data} must contain \code{ID}, \code{Time.name}, and, if \code{crwData} is hierarchical, \code{level} columns. 
#' 
#' Only rows of \code{data} with \code{ID}, \code{Time.name}, and, if \code{crwData} is hierarchical, \code{level} values that exactly match \code{crwData$crwPredict} are merged. Typically, the \code{Time.name} column
#' in \code{data} should match predicted times of locations in \code{crwData$crwPredict} (i.e. those corresponding to \code{crwData$crwPredict$locType=="p"})
#' 
#' @param crwData A \code{\link{crwData}} or \code{\link{crwHierData}} object
#' @param data A data frame containing required columns \code{ID}, \code{Time.name}, and, if \code{crwData} is hierarchical, \code{level}, plus any additional data streams and/or covariates to merge with \code{crwData}.
#' @param Time.name Character string indicating name of the time column to be used for merging
#' 
#' @return A \code{\link{crwData}} object
#' @examples
#' \dontrun{
#' # extract simulated obsData from example data
#' obsData <- miExample$obsData
#' 
#' # error ellipse model
#' err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
#'
#' # Fit crwMLE models to obsData and predict locations 
#' # at default intervals for both individuals
#' crwOut <- crawlWrap(obsData=obsData,
#'          theta=c(4,0),fixPar=c(1,1,NA,NA),
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
#' 
#' @export
crawlMerge <- function (crwData, data, Time.name) {
  UseMethod("crawlMerge")
}

#' @method crawlMerge crwData
#' @export
crawlMerge.crwData<-function(crwData, data, Time.name){
  crwFits <- crwData$crwFits
  if(!(Time.name %in% names(crwData$crwPredict))) stop(Time.name," not found in crwData$crwPredict")
  if(!("ID" %in% names(data))) stop("ID not found in data")
  if(!(Time.name %in% names(data))) stop(Time.name," not found in data")
  dataNames <- c("ID",Time.name)
  covNames <- names(data)[!(names(data) %in% dataNames)]
  #if(any(covNames %in% names(crwData$crwPredict)[!(names(crwData$crwPredict) %in% dataNames)])) stop("Data stream and/or covariate names in data cannot match column names in crwData$crwPredict")
  for(i in unique(crwData$crwPredict$ID)){
    crInd <- which(crwData$crwPredict$ID==i & crwData$crwPredict$locType=="p")
    tmp<-merge(crwData$crwPredict[crInd,!(names(crwData$crwPredict) %in% covNames)],data[which(data$ID==i),],by=Time.name,all.x=TRUE,all.y=FALSE,sort=TRUE)
    if(nrow(tmp)!=length(crInd)) stop("multiple ",Time.name," matches found for individual ",i)
    crwData$crwPredict[crInd,covNames]<-tmp[covNames]
  }
  return(crwData)
}

#' @method crawlMerge crwHierData
#' @export
crawlMerge.crwHierData<-function(crwData, data, Time.name){
  crwFits <- crwData$crwFits
  if(!(Time.name %in% names(crwData$crwPredict))) stop(Time.name," not found in crwData$crwPredict")
  if(!("ID" %in% names(data))) stop("ID not found in data")
  if(!(Time.name %in% names(data))) stop(Time.name," not found in data")
  if(!("level" %in% names(data))) stop("level not found in data")
  dataNames <- c("ID",Time.name,"level")
  covNames <- names(data)[!(names(data) %in% dataNames)]
  #if(any(covNames %in% names(crwData$crwPredict)[!(names(crwData$crwPredict) %in% dataNames)])) stop("Data stream and/or covariate names in data cannot match column names in crwData$crwPredict")
  for(i in unique(crwData$crwPredict$ID)){
    crInd <- which(crwData$crwPredict$ID==i & crwData$crwPredict$locType=="p")
    tmp<-merge(crwData$crwPredict[crInd,!(names(crwData$crwPredict) %in% covNames)],data[which(data$ID==i),],by=Time.name,all.x=TRUE,all.y=FALSE,sort=TRUE)
    if(nrow(tmp)!=length(crInd)) stop("multiple ",Time.name," matches found for individual ",i)
    crwData$crwPredict[crInd,covNames]<-tmp[covNames]
  }
  return(crwData)
}
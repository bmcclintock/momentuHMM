
#' Observation error simulation tool
#' 
#' Simulates observed location data subject to temporal irregularity and/or location measurement error
#' 
#' Simulated location data that are temporally-irregular (i.e., \code{lambda>0}) and/or with location measurement error (i.e., \code{errorEllipse!=NULL}) are returned
#' as a data frame suitable for analysis using \code{\link{crawlWrap}}.
#' 
#' @param data A \code{\link{momentuHMMData}} object with necessary field 'x' (easting/longitudinal coordinates) and 'y' (northing/latitudinal coordinates)
#' @param lambda Observation rate for location data. If \code{NULL}, location data are kept at temporally-regular intervals. Otherwise 
#' \code{lambda} is the rate parameter of the exponential distribution for the waiting times between successive location observations, i.e., 
#' \code{1/lambda} is the expected time between successive location observations. Only the 'step' and 'angle' data streams are subject to temporal irregularity;
#' any other data streams are kept at temporally-regular intervals.  Ignored unless a valid distribution for the 'step' data stream is specified.
#' @param errorEllipse List providing the bounds for the semi-major axis (\code{M}; on scale of x- and y-coordinates), semi-minor axis (\code{m}; 
#' on scale of x- and y-coordinates), and orientation (\code{r}; in degrees) of location error ellipses. If \code{NULL}, no location 
#' measurement error is simulated. If \code{errorEllipse} is specified, then each observed location is subject to bivariate normal errors as described 
#' in McClintock et al. (2015), where the components of the error ellipse for each location are randomly drawn from \code{runif(1,min(errorEllipse$M),max(errorEllipse$M))}, 
#' \code{runif(1,min(errorEllipse$m),max(errorEllipse$m))}, and \code{runif(1,min(errorEllipse$r),max(errorEllipse$r))}. If only a single value is provided for any of the 
#' error ellipse elements, then the corresponding component is fixed to this value for each location. Only the 'step' and 'angle' data streams are subject to location measurement error;
#' any other data streams are observed without error.  Ignored unless a valid distribution for the 'step' data stream is specified.
#' 
#' @return A dataframe of:
#' \item{time}{Numeric time of each observed (and missing) observation}
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{x}{Either easting or longitude observed location}
#' \item{y}{Either norting or latitude observed location}
#' \item{...}{Data streams that are not derived from location (if applicable)}
#' \item{...}{Covariates at temporally-regular true (\code{mux},\code{muy}) locations (if any)}
#' \item{mux}{Either easting or longitude true location}
#' \item{muy}{Either norting or latitude true location}
#' \item{error_semimajor_axis}{error ellipse semi-major axis (if applicable)}
#' \item{error_semiminor_axis}{error ellipse semi-minor axis (if applicable)}
#' \item{error_ellipse_orientation}{error ellipse orientation (if applicable)}
#' \item{ln.sd.x}{log of the square root of the x-variance of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' \item{ln.sd.y}{log of the square root of the y-variance of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' \item{error.corr}{correlation term of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' 
#' @seealso \code{\link{crawlWrap}}, \code{\link{prepData}}, \code{\link{simData}}
#' 
#' @references
#' McClintock BT, London JM, Cameron MF, Boveng PL. 2015. Modelling animal movement using the Argos satellite telemetry location error ellipse. 
#' Methods in Ecology and Evolution 6(3):266-277.
#' 
#' @examples 
#' # extract momentuHMMData example
#' data <- example$m$data
#' lambda <- 2 # expect 2 observations per time step
#' errorEllipse <- list(M=c(0,50),m=c(0,50),r=c(0,180))
#' obsData1 <- simObsData(data,lambda=lambda,errorEllipse=errorEllipse)
#' 
#' errorEllipse <- list(M=50,m=50,r=180)
#' obsData2 <- simObsData(data,lambda=lambda,errorEllipse=errorEllipse)
#' 
#' @importFrom mvtnorm rmvnorm
#' @importFrom argosfilter radian
#' @importFrom crawl argosDiag2Cov
#' @importFrom stats rexp
#' @export
simObsData<-function(data,lambda,errorEllipse){
  
  if(!is.momentuHMMData(data)) stop("data must be a momentuHMMData object")
  
  if(!all(c("x","y") %in% names(data)) | (is.null(lambda) & is.null(errorEllipse))){
    out <- data
  } else {
    
    if(!is.null(errorEllipse)){
      if(!is.list(errorEllipse) | any(!(c("M","m","r") %in% names(errorEllipse)))) stop("errorEllipse must be a list of scalars named 'M', 'm', and 'r'.")
      if(any(unlist(lapply(errorEllipse,length))>2)) stop('errorEllipse elements must be of length 1 or 2')
      if(any(errorEllipse$M<0 | errorEllipse$m<0)) stop("errorEllipse$M and errorEllipse$m must be >=0")
      if(any(errorEllipse$r < 0) | any(errorEllipse$r > 180)) stop('errorEllipse$r must be in [0,180]')
      
      M<-c(min(errorEllipse$M),max(errorEllipse$M))
      m<-c(min(errorEllipse$m),max(errorEllipse$m))
      r<-c(min(errorEllipse$r),max(errorEllipse$r))
      
    } else {
      M<-m<-r<-c(0,0)
    }
    if(!is.null(lambda))
      if(lambda<=0) stop('lambda must be >0')
  
    distnames<-names(data)[which(!(names(data) %in% c("ID","x","y","step","angle")))]
    
    #Get observed data based on sampling rate (lambda) and measurement error (M,m, and r)
    obsData<-data.frame()
    for(i in unique(data$ID)){
      idat<-data[which(data$ID==i),]
      nbObs<-nrow(idat)
      X<-idat$x
      Y<-idat$y
      if(!is.null(lambda)){
        t<-cumsum(c(1,rexp(2*nbObs*lambda,lambda)))
        t<-t[which(t<nbObs)]
        tind<-floor(t)
        mux<-X[tind]*(1-(t-tind))+X[tind+1]*(t-tind)
        muy<-Y[tind]*(1-(t-tind))+Y[tind+1]*(t-tind)
        mux<-c(mux,X[nbObs])
        muy<-c(muy,Y[nbObs])
        t<-c(t,nbObs)
        tind<-c(tind,nbObs)
      } else {
        t<- 1:(nbObs-1)
        tind<- floor(t)
        mux<-c(X[tind]*(1-(t-tind))+X[tind+1]*(t-tind),X[length(X)])
        muy<-c(Y[tind]*(1-(t-tind))+Y[tind+1]*(t-tind),Y[length(Y)])
        t<-c(t,nbObs)
        tind<-c(tind,nbObs)
      }
      nobs<-length(mux)
      muxy<-cbind(mux,muy)
      
      error_semimajor_axis<-tmpM<-runif(nobs,M[1],M[2])
      error_semiminor_axis<-tmpm<-runif(nobs,m[1],m[2])
      error_semimajor_axis[which(tmpM<tmpm)]<-tmpm[which(tmpM<tmpm)]
      error_semiminor_axis[which(tmpM<tmpm)]<-tmpM[which(tmpM<tmpm)]
      error_ellipse_orientation<-runif(nobs,r[1],r[2])
      rad<-argosfilter::radian(error_ellipse_orientation)
        
      #calculate bivariate normal error variance-covariance matrix (Sigma)
      sigma2x<-(error_semimajor_axis/sqrt(2))^2*sin(rad)^2+(error_semiminor_axis/sqrt(2))^2*cos(rad)^2 # x measurement error term
      sigma2y<-(error_semimajor_axis/sqrt(2))^2*cos(rad)^2+(error_semiminor_axis/sqrt(2))^2*sin(rad)^2 # y mearurement error term
      sigmaxy<-((error_semimajor_axis^2-error_semiminor_axis^2)/2)*cos(rad)*sin(rad) # measurement error ellipse rotation term
      
      xy<-t(apply(cbind(muxy,sigma2x,sigmaxy,sigma2y),1,function(x) mvtnorm::rmvnorm(1,c(x[1],x[2]),matrix(c(x[3],x[4],x[4],x[5]),2,2))))
      
      if(!is.null(errorEllipse))
        tmpobsData<-data.frame(time=t,ID=rep(i,nobs),x=xy[,1],y=xy[,2],error_semimajor_axis=error_semimajor_axis,error_semiminor_axis=error_semiminor_axis,error_ellipse_orientation=error_ellipse_orientation,crawl::argosDiag2Cov(error_semimajor_axis,error_semiminor_axis,error_ellipse_orientation),mux=mux,muy=muy)
      else
        tmpobsData<-data.frame(time=t,ID=rep(i,nobs),x=xy[,1],y=xy[,2],mux=mux,muy=muy)
      
      if(!is.null(lambda)) {
        tmpobsData<-merge(tmpobsData,data.frame(idat[,c("ID",distnames),drop=FALSE],time=1:nbObs,mux=idat$x,muy=idat$y),all = TRUE,by=c("ID","time","mux","muy"))
      } else if(length(distnames)){
        tmpobsData<-cbind(tmpobsData,idat[,distnames,drop=FALSE])
      }
      tmpobsData[order(tmpobsData$time),]
      obsData<-rbind(obsData,tmpobsData)
    }
    dnames <- names(obsData)
    out <- obsData[,c("time","ID","x","y",distnames,"mux","muy",dnames[!(dnames %in% c("time","ID","x","y",distnames,"mux","muy"))])]
  }
  return(out)
}
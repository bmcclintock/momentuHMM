#' @importFrom mvtnorm rmvnorm
#' @importFrom argosfilter radian
#' @importFrom crawl expandPred
simObsData<-function(data,dist,lambda,errorEllipse){
  
  distnames<-names(data)[which(!(names(data) %in% c("ID","x","y","step","angle")))]
  
  if(is.null(errorEllipse)) M<-m<-r<-0
  else {
    M<-errorEllipse$M
    m<-errorEllipse$m
    r<-errorEllipse$r
  }
  
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
    
    M<-runif(nobs,0,M)
    m<-runif(nobs,0,m)
    r<-runif(nobs,0,r)
    rad<-argosfilter::radian(r)
      
    #calculate bivariate normal error variance-covariance matrix (Sigma)
    sigma2x<-(M/sqrt(2))^2*sin(rad)^2+(m/sqrt(2))^2*cos(rad)^2 # x measurement error term
    sigma2y<-(M/sqrt(2))^2*cos(rad)^2+(m/sqrt(2))^2*sin(rad)^2 # y mearurement error term
    sigmaxy<-((M^2-m^2)/2)*cos(rad)*sin(rad) # measurement error ellipse rotation term
    
    #Sigma<-matrix(c(sigma2x,sigmaxy,sigmaxy,sigma2y),2,2)
    xy<-t(apply(cbind(muxy,sigma2x,sigmaxy,sigma2y),1,function(x) mvtnorm::rmvnorm(1,c(x[1],x[2]),matrix(c(x[3],x[4],x[4],x[5]),2,2))))
    
    if(!is.null(errorEllipse))
      tmpobsData<-data.frame(time=t,ID=rep(i,nobs),x=xy[,1],y=xy[,2],error_semimajor_axis=M,error_semiminor_axis=m,error_ellipse_orientation=r,mux=mux,muy=muy)
    else
      tmpobsData<-data.frame(time=t,ID=rep(i,nobs),x=xy[,1],y=xy[,2],mux=mux,muy=muy)
    
    if(!is.null(lambda)) {
      tmpobsData<-crawl::expandPred(tmpobsData,Time="time",predTime=1:nbObs,time.col=TRUE)
      tmpobsData[match(2:(nbObs-1),tmpobsData$time),c("x","y")]<-NA
      if(!is.null(errorEllipse)) tmpobsData[match(2:(nbObs-1),tmpobsData$time),c("error_semimajor_axis","error_semiminor_axis","error_ellipse_orientation")]<-NA
      tmpobsData[match(1:nbObs,tmpobsData$time),c("mux","muy")]<-cbind(X,Y)
    }
    if(length(distnames)){
      addStreams<-matrix(NA,nrow=nrow(tmpobsData),ncol=length(distnames))
      colnames(addStreams)<-distnames
      tmpobsData<-cbind(tmpobsData,addStreams)
      tmpobsData[match(1:nbObs,tmpobsData$time),distnames]<-as.data.frame(idat)[,distnames,drop=FALSE]
    }
    obsData<-rbind(obsData,tmpobsData)
  }
  obsData
}
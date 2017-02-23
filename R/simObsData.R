#' @importFrom mvtnorm rmvnorm
#' @importFrom argosfilter radian
#' @importFrom crawl expandPred argosDiag2Cov
#' @importFrom stats rexp
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
    
    error_semimajor_axis<-tmpM<-runif(nobs,0,M)
    error_semiminor_axis<-tmpm<-runif(nobs,0,m)
    error_semimajor_axis[which(tmpM<tmpm)]<-tmpm[which(tmpM<tmpm)]
    error_semiminor_axis[which(tmpM<tmpm)]<-tmpM[which(tmpM<tmpm)]
    error_ellipse_orientation<-runif(nobs,0,r)
    rad<-argosfilter::radian(error_ellipse_orientation)
      
    #calculate bivariate normal error variance-covariance matrix (Sigma)
    sigma2x<-(error_semimajor_axis/sqrt(2))^2*sin(rad)^2+(error_semiminor_axis/sqrt(2))^2*cos(rad)^2 # x measurement error term
    sigma2y<-(error_semimajor_axis/sqrt(2))^2*cos(rad)^2+(error_semiminor_axis/sqrt(2))^2*sin(rad)^2 # y mearurement error term
    sigmaxy<-((error_semimajor_axis^2-error_semiminor_axis^2)/2)*cos(rad)*sin(rad) # measurement error ellipse rotation term
    
    #Sigma<-matrix(c(sigma2x,sigmaxy,sigmaxy,sigma2y),2,2)
    xy<-t(apply(cbind(muxy,sigma2x,sigmaxy,sigma2y),1,function(x) mvtnorm::rmvnorm(1,c(x[1],x[2]),matrix(c(x[3],x[4],x[4],x[5]),2,2))))
    
    if(!is.null(errorEllipse))
      tmpobsData<-data.frame(time=t,ID=rep(i,nobs),x=xy[,1],y=xy[,2],error_semimajor_axis=error_semimajor_axis,error_semiminor_axis=error_semiminor_axis,error_ellipse_orientation=error_ellipse_orientation,crawl::argosDiag2Cov(error_semimajor_axis,error_semiminor_axis,error_ellipse_orientation),mux=mux,muy=muy)
    else
      tmpobsData<-data.frame(time=t,ID=rep(i,nobs),x=xy[,1],y=xy[,2],mux=mux,muy=muy)
    
    if(!is.null(lambda)) {
      tmpobsData<-crawl::expandPred(tmpobsData,Time="time",predTime=1:nbObs,time.col=TRUE)
      tmpobsData[match(2:(nbObs-1),tmpobsData$time),c("x","y")]<-NA
      if(!is.null(errorEllipse)) tmpobsData[match(2:(nbObs-1),tmpobsData$time),c("error_semimajor_axis","error_semiminor_axis","error_ellipse_orientation","ln.sd.x","ln.sd.y","error.corr","diag.check")]<-NA
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
  dnames<-names(obsData)
  obsData[,c("time","ID","x","y",distnames,"mux","muy",dnames[!(dnames %in% c("time","ID","x","y",distnames,"mux","muy"))])]
}
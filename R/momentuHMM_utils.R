momentuHMMdists<-c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm')
angledists<-c('wrpcauchy','vm')
stepdists<-c('gamma','weibull','exp','lnorm')
singleParmdists<-c('exp','pois')
nonnegativedists<-c('gamma','weibull','exp','lnorm','pois')
zeroInflationdists<-c('gamma','weibull','exp','lnorm','beta')

getboundInd<-function(DM=NULL){
  if(is.null(DM)){
    Ind<-NULL
  } else {
    Ind<-apply(apply(DM,1,function(x) apply(unique(DM),1,function(y) all(y==x))),2,function(x) which(x))
  }
  Ind
}

n2wDM<-function(bounds,DM,par,cons,workcons){

  a<-bounds[,1]
  b<-bounds[,2]

  piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind1<-which(piInd)
  ind2<-which(!piInd)
  
  p<-numeric(nrow(DM))
  
  if(length(ind1)) p[ind1] <- (tan(par[ind1]/2)-workcons[ind1])^(1/cons[ind1])
  
  p[ind2] <- par[ind2]
  
  ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
  ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
  ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
  
  p[ind21]<-(log(par[ind21]-a[ind21])-workcons[ind21])^(1/cons[ind21])
  p[ind22]<-(logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))-workcons[ind22])^(1/cons[ind22])
  p[ind23]<-(-log(-par[ind23]+b[ind23])-workcons[ind23])^(1/cons[ind23])
  
  p
}     

w2nDM<-function(wpar,bounds,DM,DMind,cons,workcons,nbObs,k=0){
  
  Par<-numeric(length(wpar))

  a<-bounds[,1]
  b<-bounds[,2]

  piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind1<-which(piInd)
  ind2<-which(!piInd)
  
  XB <- p <- getXB(DM,nbObs,wpar,cons,workcons,DMind)
 
  if(length(ind1))
    p[ind1,] <- (2*atan(XB[ind1,]))
  
  ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
  ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
  ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
  
  p[ind21,] <- (exp(XB[ind21,,drop=FALSE])+a[ind21])
  p[ind22,] <- ((b[ind22]-a[ind22])*boot::inv.logit(XB[ind22,,drop=FALSE])+a[ind22])
  p[ind23,] <- -(exp(-XB[ind23,,drop=FALSE]) - b[ind23])
  
  if(any(p<a | p>b))
    stop("Scaling error. Check initial values and bounds.")
  
  if(k) {
    p <- p[k]
    return(p)
  } else {
    if(DMind) p <- matrix(p,length(ind1)+length(ind2),nbObs)
    return(list(p=p,Par=Par))
  }
}   
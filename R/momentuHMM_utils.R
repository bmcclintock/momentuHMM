getboundInd<-function(DM=NULL){
  if(is.null(DM)){
    Ind<-NULL
  } else {
    Ind<-apply(apply(DM,1,function(x) apply(unique(DM),1,function(y) all(y==x))),2,function(x) which(x))
  }
  Ind
}

n2wDM<-function(bounds,DM,par,cons,logitcons){

  a<-bounds[,1]
  b<-bounds[,2]

  piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind1<-which(piInd)
  ind2<-which(!piInd)
  
  p<-numeric(nrow(DM))
  
  if(length(ind1)) p[ind1] <- (solve(unique(DM)[ind1,ind1],tan(par[ind1]/2))-logitcons[ind1])^(1/cons[ind1])
  par<-par[ind2]
  a<-a[ind2]
  b<-b[ind2]
  cons<-cons[ind2]
  logitcons<-logitcons[ind2]
  DM<-DM[ind2,ind2]
  if(all(is.finite(a)) & all(is.infinite(b))){
    p[ind2]<-solve(unique(DM),log(par-a))^(1/cons)
  } else if(all(is.finite(a)) & all(is.finite(b))){
    p[ind2]<-(solve(unique(DM),logit((par-a)/(b-a)))-logitcons)^(1/cons)
  } else {
    ind21<-which(is.infinite(b))
    ind22<-which(is.finite(b))
    p[ind2][ind21]<-solve(unique(DM)[ind21,ind21],log(par[ind21]-a[ind21]))^(1/cons[ind21])
    p[ind2][ind22]<-(solve(unique(DM)[ind22,ind22],logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))))^(1/cons[ind22])
  }
  p
}     

w2nDM<-function(wpar,bounds,DM,DMind,cons,logitcons,nbObs,parSize,k=0){
  
  Par<-numeric(length(wpar))
  
  #ord<-order(rep(1:(nrow(DM)/parSize),parSize))
  #bounds<-bounds[ord,]

  a<-bounds[,1]
  b<-bounds[,2]

  piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind1<-which(piInd)
  ind2<-which(!piInd)
  
  XB <- p <- getXB(DM,nbObs,wpar,cons,logitcons,DMind)
  
  #p <- matrix(0,length(ind1)+length(ind2),nbObs)
 
  if(length(ind1))
    p[ind1,] <- (2*atan(XB[ind1,]))
  
  a<-a[ind2]
  b<-b[ind2]
  
  if(all(is.finite(a[ind2])) & all(is.infinite(b[ind2]))){
    p[ind2,]<-(exp(XB[ind2,,drop=FALSE])+a[ind2])
  } else if(all(is.finite(a[ind2])) & all(is.finite(b[ind2]))){
    p[ind2,] <- ((b[ind2]-a[ind2])*boot::inv.logit(XB[ind2,,drop=FALSE])+a[ind2])
  } else {
    ind21<-which(is.infinite(b[ind2]))
    ind22<-which(is.finite(b[ind2]))
    p[ind2[ind21],]<-(exp(XB[ind2[ind21],,drop=FALSE])+a[ind2[ind21]])
    p[ind2[ind22],]<-((b[ind2[ind22]]-a[ind2[ind22]])*boot::inv.logit(XB[ind2[ind22],,drop=FALSE])+a[ind2[ind22]])
  }
  
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
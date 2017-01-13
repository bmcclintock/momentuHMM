getXB<-function(DM,nbObs,wpar,cons,logitcons,DMind){
  Xvec<-wpar^cons+logitcons
  if(DMind){
    XB <- matrix(Xvec%*%t(DM),nrow(DM),1)
  } else {
    nr<-nrow(DM)
    nc<-ncol(DM)
    XB<-matrix(0,nrow(DM),nbObs)
    for(i in 1:nr){
      DMrow<-DM[i,]
      for(j in 1:nc){
        XB[i,]<-XB[i,]+DMrow[[j]]*Xvec[j]
      }
    }
  }
  XB
}
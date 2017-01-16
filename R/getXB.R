getXB<-function(DM,nbObs,wpar,cons,workcons,DMind){
  Xvec<-wpar^cons+workcons
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
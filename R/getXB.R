getXB<-function(DM,nbObs,wpar,cons,workcons,DMind){
  Xvec<-wpar^cons
  if(DMind){
    XB <- matrix(Xvec%*%t(DM)+workcons,nrow(DM),1)
  } else {
    nr<-nrow(DM)
    nc<-ncol(DM)
    XB<-matrix(workcons,nrow(DM),nbObs)
    for(i in 1:nr){
      DMrow<-DM[i,]
      for(j in 1:nc){
        XB[i,]<-XB[i,]+DMrow[[j]]*Xvec[j]
      }
    }
  }
  XB
}
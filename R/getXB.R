getXB<-function(DM,nbObs,wpar,cons,logitcons,ord,DMind){
  Xvec<-wpar^cons+logitcons
  if(DMind){
    XB <- matrix(((wpar^cons+logitcons)%*%t(DM))[ord],nrow(DM),1)
  } else {
    nr<-nrow(DM)
    nc<-ncol(DM)
    XB<-matrix(0,nrow(DM),nbObs)
    for(i in 1:nr){
      DMrow<-DM[i,]
      for(j in 1:nc){
        XB[ord[i],]<-XB[ord[i],]+DMrow[[j]]*Xvec[j]
      }
    }
  }
  XB
}
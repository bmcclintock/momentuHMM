getXB<-function(DM,nbObs,wpar,cons,logitcons){
  Xvec<-wpar^cons+logitcons
  nr<-nrow(DM)
  nc<-ncol(DM)
  XB<-matrix(0,nrow(DM),nbObs)
  for(i in 1:nr){
    DMrow<-DM[i,]
    for(j in 1:nc){
      XB[i,]<-XB[i,]+DMrow[[j]]*Xvec[j]
    }
  }
  XB
}
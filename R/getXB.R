getXB<-function(DM,nbObs,wpar,cons,logitcons,ord){
  Xvec<-wpar^cons+logitcons
  nr<-nrow(DM)
  nc<-ncol(DM)
  XB<-matrix(0,nrow(DM),nbObs)
  for(i in 1:nr){
    DMrow<-DM[i,]
    for(j in 1:nc){
      XB[ord[i],]<-XB[ord[i],]+DMrow[[j]]*Xvec[j]
    }
  }
  XB
}
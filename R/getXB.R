getXB<-function(DM,nbObs,wpar,cons,workcons,DMind,circularAngleMean,nbStates){
  Xvec<-wpar^cons+workcons
  if(!circularAngleMean){
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
  } else {
    #if(DMind){
    #  meanind<-unique(unlist(apply(DM[1:nbStates,,drop=FALSE],1,function(x) which(x!=0))))
    #  XB <- matrix(0,nrow(DM),1)
    #  XB[meanind,] <- atan2(Xvec%*%t(sin(DM[meanind,])),1+Xvec%*%t(cos(DM[meanind,])))[meanind]
    #  XB[nbStates+1:nbStates,] <- matrix((Xvec%*%t(DM))[nbStates+1:nbStates],nbStates,1)
    #} else {
      meanind<-which((apply(DM[1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      nc<-apply(DM,1:2,function(x) !all(unlist(x)==0))
      XB<-matrix(0,nrow(DM),nbObs)
      XB1<-XB2<-matrix(0,nbStates,nbObs)
      for(i in meanind){
        DMrow<-DM[i,]
        for(j in which(nc[i,])){
          XB1[i,]<-XB1[i,]+sin(DMrow[[j]])*Xvec[j]
          XB2[i,]<-XB2[i,]+cos(DMrow[[j]])*Xvec[j]
        }
        XB[i,] <- atan2(XB1[i,],1+XB2[i,])
      }
      for(i in nbStates+1:nbStates){
        DMrow<-DM[i,]
        for(j in which(nc[i,])){
          XB[i,]<-XB[i,]+DMrow[[j]]*Xvec[j]
        }
      }
    #}
  }
  XB
}
getXB<-function(DM,nbObs,wpar,cons,workcons,DMind,circularAngleMean,consensus,nbStates,nc,meanind){
  Xvec<-wpar^cons+workcons
  nrw <- nrow(DM)
  ncl <- ncol(DM)
  if(isFALSE(circularAngleMean)){
    if(DMind){
      XB <- DM%*%Xvec
    } else {
      XB <- XBloop_rcpp(DM, Xvec, nbObs, nrw, ncl, FALSE, FALSE, 1:nrw-1, nc, nbStates)
    }
  } else {
    if(length(meanind)) {
      XB<-XBloop_rcpp(DM,Xvec,nbObs,nrw,ncl,TRUE,consensus,meanind-1,nc,nbStates,as.numeric(circularAngleMean))
      if(consensus) {
        l_t <- matrix(1,nrow(XB),ncol(XB))
        l_t[nbStates+meanind,] <- XB[nbStates+meanind,]
      }
    } else {
      XB <- matrix(0,nrw,nbObs)
      if(consensus) l_t <- matrix(1,nrow(XB),ncol(XB))
    }
    cInd <- which(colSums(nc[nbStates+1:nbStates,,drop=FALSE])>0)
    Xvec2 <- Xvec[cInd - (ncl-length(Xvec))]
    XB[nbStates+1:nbStates,]<-XBloop_rcpp(DM[nbStates+1:nbStates,cInd,drop=FALSE],Xvec2,nbObs,nbStates,length(Xvec2),FALSE,FALSE,1:nbStates-1,nc[nbStates+1:nbStates,cInd,drop=FALSE],nbStates)
  }
  if(consensus){
    return(list(XB=XB,l_t=l_t))
  } else return(XB)
}
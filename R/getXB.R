getXB<-function(DM,nbObs,wpar,cons,workcons,DMind,circularAngleMean,nbStates,nc,meanind){
  Xvec<-wpar^cons+workcons
  nrw <- nrow(DM)
  ncl <- ncol(DM)
  if(!circularAngleMean){
    if(DMind){
      XB <- DM%*%Xvec
    } else {
      XB <- XBloop_rcpp(DM, Xvec, nbObs, nrw, ncl, circularAngleMean, 1:nrw-1, nc)
    }
  } else {
    if(length(meanind)) XB<-XBloop_rcpp(DM,Xvec,nbObs,nrw,ncl,circularAngleMean,meanind-1,nc)
    else XB <- matrix(0,nrw,nbObs)
    XB[nbStates+1:nbStates,]<-XBloop_rcpp(DM[nbStates+1:nbStates,],Xvec,nbObs,nbStates,ncl,FALSE,1:nbStates-1,nc[nbStates+1:nbStates,])
  }
  XB
}
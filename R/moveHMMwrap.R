# wrapper for calling moveHMM::fitHMM when appropriate models are specified
#' @importFrom moveHMM fitHMM
moveHMMwrap<-function(data,nbStates,dist,Par,beta0,delta0,estAngleMean,formula,stationary,nlmPar,fit,nbAnimals){
  data <- moveData(data)
  distnames<-names(dist)
  
  verbose <- ifelse(is.null(nlmPar$print.level),0,nlmPar$print.level)
  nlmPar$print.level <- NULL
  
  if(any(distnames %in% "angle")){
    
    angleMean<-rep(0,nbStates)
    if(!is.null(estAngleMean$angle)) 
      if(estAngleMean$angle) angleMean<-NULL
      
    out<-moveHMM::fitHMM(data,nbStates,Par$step,Par$angle,beta0,delta0,formula,dist$step,dist$angle,angleMean,stationary,knownStates=NULL,verbose,nlmPar,fit)
    estAngleMean<-list(step=FALSE,angle=out$conditions$estAngleMean)
    zeroInflation<-list(step=out$conditions$zeroInflation,angle=FALSE)
    cons<-list(step=rep(1,length(Par$step)),angle=rep(1,length(Par$angle)))
    workcons<-list(step=rep(0,length(Par$step)),angle=rep(0,length(Par$angle)))
    mle<-out$mle
    mle$angle<-matrix(c(t(mle$anglePar)),length(mle$anglePar),1)
  } else {
    out<-moveHMM::fitHMM(data,nbStates,Par$step,Par$angle,beta0,delta0,formula,dist$step,"none",NULL,stationary,knownStates=NULL,verbose,nlmPar,fit)
    estAngleMean<-list(step=FALSE)
    zeroInflation<-list(step=out$conditions$zeroInflation)
    cons<-list(step=rep(1,length(Par$step)))
    workcons<-list(step=rep(0,length(Par$step)))
    mle<-out$mle
  }
  mod<-out$mod
  mle$step<-matrix(c(t(mle$stepPar)),length(mle$stepPar),1)
  mle$delta<-matrix(mle$delta,nrow=nbAnimals,ncol=nbStates,byrow=TRUE)
  return(list(mod=mod,mle=mle))
}
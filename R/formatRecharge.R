formatRecharge <- function(m,data,covs=NULL,par=m$mle){
  
  nbStates <- length(m$stateNames)
  
  formula <- m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  aInd <- NULL
  nbAnimals <- length(unique(data$ID))
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    
    g0covs <- model.matrix(recharge$g0,data[aInd,])
    nbG0covs <- ncol(g0covs)-1
    recovs <- model.matrix(recharge$theta,data)
    nbRecovs <- ncol(recovs)-1
    data$recharge<-rep(0,nrow(data))
    for(i in 1:nbAnimals){
      idInd <- which(data$ID==unique(data$ID)[i])
      if(nbRecovs){
        g0 <- par$g0 %*% t(g0covs[i,,drop=FALSE])
        theta <- par$theta
        data$recharge[idInd] <- cumsum(c(g0,theta%*%t(recovs[idInd[-length(idInd)],])))
      }
    }
    newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
    if(!is.null(covs)){
      if(is.null(covs$recharge)) covs$recharge <- mean(data$recharge)
    }
  } else {
    nbG0covs <- 0
    nbRecovs <- 0
    g0covs <- NULL
    recovs <- NULL
  }
  
  tmpCovs <- model.matrix(newformula,data)
  if(is.null(covs)) {
    covs <- tmpCovs
  }
  nbCovs <- ncol(tmpCovs)-1 # substract intercept column
  
  list(data=data,newformula=newformula,recharge=recharge,covs=covs,nbCovs=nbCovs,nbG0covs=nbG0covs,nbRecovs=nbRecovs,g0covs=g0covs,recovs=recovs,aInd=aInd)
}
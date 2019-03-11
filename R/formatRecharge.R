formatRecharge <- function(nbStates,formula,data,covs=NULL,par=NULL,hierRecharge=NULL){
  
  nbAnimals <- length(unique(data$ID))
  
  if(is.null(hierRecharge)){
    newForm <- newFormulas(formula,nbStates)
    newformula <- newForm$newformula
    recharge <- newForm$recharge
    rechargeName <- "recharge"
    
    aInd <- NULL
    for(i in 1:nbAnimals){
      aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
    }
    
    if(!is.null(recharge)){
      g0covs <- model.matrix(recharge$g0,data[aInd,])
      recovs <- model.matrix(recharge$theta,data)
    }
  } else {
    newformula <- formula
    recharge <- hierRecharge[[1]]
    #rechargeName <- paste0("recharge",gsub("level","",names(hierRecharge)))
    rechargeName <- "recharge"
    
    aInd <- NULL
    for(i in 1:nbAnimals){
      aInd <- c(aInd,which(data$ID==unique(data$ID)[i] & data$level==gsub("level","",names(hierRecharge)))[1])
    }
    
    if(!is.null(recharge)){
      g0covs <- model.matrix(recharge$g0,data[aInd,])
      recovs <- model.matrix(recharge$theta,data)
      recovs[which(data$level!=gsub("level","",names(hierRecharge))),] <- 0
    }
  }
  
  if(!is.null(recharge)){
    
    nbG0covs <- ncol(g0covs)-1
    nbRecovs <- ncol(recovs)-1
    data[[rechargeName]]<-rep(0,nrow(data))
    
    if(is.null(par)) par <- list(g0=rep(0,nbG0covs+1),theta=rep(0,nbRecovs+1))
    for(i in 1:nbAnimals){
      idInd <- which(data$ID==unique(data$ID)[i])
      if(nbRecovs){
        g0 <- par$g0 %*% t(g0covs[i,,drop=FALSE])
        theta <- par$theta
        data[[rechargeName]][idInd] <- cumsum(c(g0,theta%*%t(recovs[idInd[-length(idInd)],])))
      }
    }
    if(is.null(hierRecharge)) newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
    if(!is.null(covs)){
      if(is.null(covs[,colnames(covs)[which(grepl(rechargeName,colnames(covs)))]])) covs[,colnames(covs)[which(grepl(rechargeName,colnames(covs)))]] <- mean(data[[rechargeName]])
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
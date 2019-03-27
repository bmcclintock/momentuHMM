formatRecharge <- function(nbStates,formula,data,covs=NULL,par=NULL){
  
  nbAnimals <- length(unique(data$ID))
  dataNames <- colnames(data)
  
  newForm <- newFormulas(formula,nbStates)
  newformula <- newForm$newformula
  formulaStates <- newForm$formulaStates
  recharge <- newForm$recharge
  rechargeName <- "recharge"
  
  aInd <- NULL
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    if(inherits(data,"hierarchical")){
      tmpg0covs <- model.matrix(recharge$g0,data)
      g0covs <- tmpg0covs[aInd,,drop=FALSE]
      for(iLevel in levels(data$level)){
        lInd <- paste0("I((level == \"",iLevel,"\") * 1)")
        if(any(lInd %in% rownames(attr(terms(recharge$g0),"factors")))){
          for(i in 1:nbAnimals){
            #reFact <- attr(terms(recharge$g0),"factors")[lInd,,drop=FALSE]
            g0covs[i,] <- tmpg0covs[which(data$level==iLevel & data$ID==unique(data$ID)[i])[1],,drop=FALSE]
          }
        }
      }
    } else g0covs <- model.matrix(recharge$g0,data[aInd,])
    recovs <- model.matrix(recharge$theta,data)
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
    #if(is.null(hierRecharge)) newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
    if(!is.null(covs)){
      covs[[rechargeName]] <- mean(data[[rechargeName]])
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
  
  list(newdata=data[,which(!(names(data) %in% dataNames)),drop=FALSE],newformula=newformula,formulaStates=formulaStates,recharge=recharge,covs=covs,nbCovs=nbCovs,nbG0covs=nbG0covs,nbRecovs=nbRecovs,g0covs=g0covs,recovs=recovs,aInd=aInd)
}
formatRecharge <- function(nbStates,formula,betaRef,data,covs=NULL,par=NULL){
  
  nbAnimals <- length(unique(data$ID))
  dataNames <- colnames(data)
  
  newForm <- newFormulas(formula,nbStates,betaRef)
  newformula <- newForm$newformula
  formulaStates <- newForm$formulaStates
  recharge <- hierRecharge <- newForm$recharge
  
  aInd <- NULL
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    
    recharge <- expandRechargeFormulas(hierRecharge)
    
    if(inherits(data,"hierarchical")){
      
      recLevels <- length(hierRecharge)
      recLevelNames <- names(hierRecharge)
      
      tmpg0covs <- model.matrix(recharge$g0,data)
      g0covs <- tmpg0covs[rep(aInd,each=recLevels)+rep(0:(recLevels-1),nbAnimals),,drop=FALSE]
      for(i in 1:nbAnimals){
        for(iLevel in 1:recLevels){
          g0covs[(i-1)*recLevels+iLevel,] <- tmpg0covs[which(data$level==iLevel & data$ID==unique(data$ID)[i])[1],,drop=FALSE]
        }
      }
      recovs <- model.matrix(recharge$theta,data)
      rechargeNames <- paste0("recharge",gsub("level","",recLevelNames))
      colInd <- lapply(recLevelNames,function(x) which(grepl(paste0("I((level == \"",gsub("level","",x),"\")"),colnames(recovs),fixed=TRUE)))
    } else {
      g0covs <- model.matrix(recharge$g0,data[aInd,])
      recovs <- model.matrix(recharge$theta,data)
      recLevels <- 1
      rechargeNames <- "recharge"
      colInd <- list(1:ncol(recovs))
    }
    
    nbG0covs <- ncol(g0covs)-1
    nbRecovs <- ncol(recovs)-1
    data[rechargeNames]<- list(rep(0,nrow(data)))
    
    if(is.null(par)) par <- list(g0=rep(0,nbG0covs+1),theta=rep(0,nbRecovs+1))
    for(i in 1:nbAnimals){
      for(iLevel in 1:recLevels){
        idInd <- which(data$ID==unique(data$ID)[i])
        if(nbRecovs){
          if(!all(names(par$theta)==colnames(recovs)) | !all(names(par$g0)==colnames(g0covs))) stop("column name mismatch in hierarchical recharge model -- please report to brett.mcclintock@noaa.gov")
          g0 <- par$g0 %*% t(g0covs[(i-1)*recLevels+iLevel,,drop=FALSE])
          theta <- par$theta
          data[[rechargeNames[iLevel]]][idInd] <- cumsum(c(g0,theta[colInd[[iLevel]]]%*%t(recovs[idInd[-length(idInd)],colInd[[iLevel]]])))
        }
      }
    }
    #if(is.null(hierRecharge)) newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
    if(!is.null(covs)){
      covs[rechargeNames] <- lapply(data[rechargeNames],mean)
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
  
  list(newdata=data[,which(!(names(data) %in% dataNames)),drop=FALSE],newformula=newformula,formulaStates=formulaStates,recharge=recharge,hierRecharge=hierRecharge,covs=covs,nbCovs=nbCovs,nbG0covs=nbG0covs,nbRecovs=nbRecovs,g0covs=g0covs,recovs=recovs,aInd=aInd)
}

expandRechargeFormulas <- function(recharge){
  
  rec <- recharge
  
  if(!is.null(recharge)){
    recLevels <- length(recharge)
    
    if(!all(names(recharge) %in% c("g0","theta"))){
      g0forms <- lapply(recharge,function(x) x$g0)
      thetaforms <- lapply(recharge,function(x) x$theta)
      
      g0form <- thetaform <- NULL
      for(i in 1:recLevels){
        g0form <- c(g0form,as.character(g0forms[[i]][-1]))
        thetaform <- c(thetaform,as.character(thetaforms[[i]][-1]))
      }
      g0form <- as.formula(paste0("~",paste0(g0form,collapse="+")))
      thetaform <- as.formula(paste0("~",paste0(thetaform,collapse="+")))
    } else {
      g0form <- recharge$g0
      thetaform <- recharge$theta
    }
    rec <- list(g0=g0form,theta=thetaform)
  }
  rec
}
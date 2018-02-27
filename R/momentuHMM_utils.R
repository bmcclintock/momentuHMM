momentuHMMdists<-sort(c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm','norm','bern',"vmConsensus"))
moveHMMdists<-sort(c('gamma','weibull','exp','lnorm','wrpcauchy','vm'))
angledists<-sort(c('wrpcauchy','vm','vmConsensus'))
stepdists<-sort(c('gamma','weibull','exp','lnorm'))
singleParmdists<-sort(c('exp','pois','bern'))
nonnegativedists<-sort(c('gamma','weibull','exp','lnorm','pois'))
zeroInflationdists<-sort(c('gamma','weibull','exp','lnorm','beta'))
oneInflationdists<-sort(c('beta'))
integerdists<-sort(c('bern','pois'))
splineList<-c("bs","ns","bSpline","mSpline","cSpline","iSpline")
meansList<-c("numeric","integer","logical","Date","POSIXlt","POSIXct","difftime")
meansListNoTime<-c("numeric","integer","logical")
plotArgs <- c("cex","cex.main","cex.lab","cex.axis","cex.legend","lwd","asp","legend.pos")
fitMethods<-c("nlm","Nelder-Mead","SANN")

# this function maintains backwards compatibility with momentuHMM versions <1.1.2 (formulaDelta) and <1.4.0 (workBounds)
delta_bc <- function(m){
  
  if(is.momentuHMM(m) | is.miSum(m)){
    if(is.null(m$conditions$formulaDelta)) {
      m$conditions$formulaDelta <- ~1
      aInd <- NULL
      nbAnimals <- length(unique(m$data$ID))
      for(i in 1:nbAnimals)
        aInd <- c(aInd,which(m$data$ID==unique(m$data$ID)[i])[1])
      m$covsDelta <- model.matrix(m$conditions$formulaDelta,m$data[aInd,,drop=FALSE]) 
      if(is.miSum(m)){
        m$Par$real$delta$est <- matrix(m$Par$real$delta$est,nrow=nrow(m$covsDelta),ncol=length(m$stateNames),byrow=TRUE,dimnames = list(paste0("ID:",unique(m$data$ID)),m$stateNames))
      } else {
        m$mle$delta <- matrix(m$mle$delta,nrow=nrow(m$covsDelta),ncol=length(m$stateNames),byrow=TRUE,dimnames = list(paste0("ID:",unique(m$data$ID)),m$stateNames))
        #m$CIreal <- CIreal(m)
        m$CIbeta <- CIbeta(m)
      }
    }
    if(is.null(m$conditions$workBounds)){
      distnames <- names(m$conditions$dist)
      
      parCount<- lapply(m$conditions$fullDM,ncol)
      for(i in distnames[unlist(m$conditions$circularAngleMean)]){
        parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
      }
      parindex <- c(0,cumsum(unlist(parCount))[-length(m$conditions$fullDM)])
      names(parindex) <- distnames
      
      workBounds <- vector('list',length(distnames))
      names(workBounds) <- distnames
      if(is.miSum(m)){
        beta <- m$Par$beta$beta$est
        delta <- m$Par$beta$delta$est
      } else {
        beta <- m$CIbeta$beta$est
        delta <- m$CIbeta$delta$est
      }
      m$conditions$workBounds <- getWorkBounds(workBounds,distnames,m$mod$estimate,parindex,parCount,m$conditions$DM,beta,delta)
    }

  } else if(!is.miHMM(m) & any(unlist(lapply(m,is.momentuHMM)))){
    m <- HMMfits(m)
  }
  m
}
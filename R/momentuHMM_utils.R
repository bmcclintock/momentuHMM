momentuHMMdists<-sort(c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm','norm','bern'))
angledists<-sort(c('wrpcauchy','vm'))
stepdists<-sort(c('gamma','weibull','exp','lnorm'))
singleParmdists<-sort(c('exp','pois','bern'))
nonnegativedists<-sort(c('gamma','weibull','exp','lnorm','pois'))
zeroInflationdists<-sort(c('gamma','weibull','exp','lnorm','beta'))
oneInflationdists<-sort(c('beta'))
integerdists<-sort(c('bern','pois'))
splineList<-c("bs","ns","bSpline","mSpline","cSpline","iSpline")
meansList<-c("numeric","integer","logical","Date","POSIXlt","POSIXct","difftime")
meansListNoTime<-c("numeric","integer","logical")

# this function maintains backwards compatibility with momentuHMM versions <1.1.2
delta_bc <- function(m){
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
  m
}
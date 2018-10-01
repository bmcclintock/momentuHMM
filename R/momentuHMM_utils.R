momentuHMMdists<-sort(c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm','norm','bern','vmConsensus','mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3','rw_norm'))
moveHMMdists<-sort(c('gamma','weibull','exp','lnorm','wrpcauchy','vm'))
angledists<-sort(c('wrpcauchy','vm','vmConsensus'))
stepdists<-sort(c('gamma','weibull','exp','lnorm'))
singleParmdists<-sort(c('exp','pois','bern'))
nonnegativedists<-sort(c('gamma','weibull','exp','lnorm','pois'))
zeroInflationdists<-sort(c('gamma','weibull','exp','lnorm','beta'))
oneInflationdists<-sort(c('beta'))
integerdists<-sort(c('bern','pois'))
mvndists <- c('mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3')
rwdists <- c('rw_norm','rw_mvnorm2','rw_mvnorm3')
splineList<-c("bs","ns","bSpline","mSpline","cSpline","iSpline")
meansList<-c("matrix","numeric","integer","logical","Date","POSIXlt","POSIXct","difftime")
meansListNoTime<-c("numeric","integer","logical")
plotArgs <- c("cex","cex.main","cex.lab","cex.axis","cex.legend","lwd","asp","legend.pos")
fitMethods<-c("nlm","Nelder-Mead","SANN")

dmvnorm2 <- dmvnorm3 <- drw_mvnorm2 <- drw_mvnorm3 <- function(x,mean,sigma){
  dmvnorm_rcpp(x,mean,sigma)
}

drw_norm <- stats::dnorm

RWdata <- function(dist,data){
  distnames <- names(dist)
  if(any(unlist(dist) %in% rwdists)){
    colInd <- NULL
    for(i in distnames){
      if(dist[[i]] %in% rwdists){
        if(dist[[i]] %in% mvndists){
          data[[paste0(i,".x_tm1")]] <- data[[paste0(i,".x")]]
          data[[paste0(i,".y_tm1")]] <- data[[paste0(i,".y")]]
          if(dist[[i]]=="rw_mvnorm2"){
            colInd <- unique(c(colInd,colnames(data)[which(!(colnames(data) %in% c(distnames[which(!(distnames %in% i))],paste0(i,".x"),paste0(i,".y"))))]))
            #data[,colInd] <- rbind(rep(NA,length(colInd)),data[-nrow(data),colInd])
          } else if(dist[[i]]=="rw_mvnorm3"){
            data[[paste0(i,".z_tm1")]] <- data[[paste0(i,".z")]]
            colInd <- unique(c(colInd,colnames(data)[which(!(colnames(data) %in% c(distnames[which(!(distnames %in% i))],paste0(i,".x"),paste0(i,".y"),paste0(i,".z"))))]))
            #data[,colInd] <- rbind(rep(NA,length(colInd)),data[-nrow(data),colInd])
          }
        } else {
          data[[paste0(i,"_tm1")]] <- data[[i]]
          colInd <- unique(c(colInd,colnames(data)[which(!(colnames(data) %in% distnames))]))
          #data[,colInd] <- rbind(rep(NA,length(colInd)),data[-nrow(data),colInd])
        }
      }
    }
    data[,colInd] <- rbind(rep(NA,length(colInd)),data[-nrow(data),colInd])
    data <- data[-1,,drop=FALSE]
  }
  data
}

# startup message
#' @importFrom utils packageDescription available.packages
print.momentuHMM.version <- function()
{ pkgDescr <- utils::packageDescription("momentuHMM")
  hello <- paste("momentuHMM ",pkgDescr$Version," (",pkgDescr$Date,")",sep="")
  curVersion <- tryCatch(suppressWarnings(utils::available.packages(repos = "http://cran.us.r-project.org")["momentuHMM","Version"]),error=function(e) e)
  packageStartupMessage(hello)
  if(!inherits(curVersion,"error")){
    if(pkgDescr$Version<curVersion) warning("  A newer version (",curVersion,") is available from CRAN")
  }
}

.onAttach <- function(...) { 
  print.momentuHMM.version()
}

# suppress RNG warning when using %dorng%
muffleRNGwarning <- function(w) {
  if(any(grepl("Foreach loop had changed the current RNG type: RNG was restored to same type, next state",w)))
    invokeRestart("muffleWarning")
}

#' @importFrom MASS ginv
# this function maintains backwards compatibility with momentuHMM versions <1.1.2 (formulaDelta), <1.4.0 (workBounds), and <1.4.3 (betaCons)
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
      for(i in distnames[!unlist(lapply(m$conditions$circularAngleMean,isFALSE))]){
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
      beta <- list(beta=beta,g0=m$mle$g0,theta=m$mle$theta)
      m$conditions$workBounds <- getWorkBounds(workBounds,distnames,m$mod$estimate,parindex,parCount,m$conditions$DM,beta,delta)
    }
    if(is.null(m$conditions$betaCons)){
      if(is.miSum(m) & !is.null(m$Par$beta$beta)) m$conditions$betaCons <- matrix(1:length(m$Par$beta$beta$est),nrow(m$Par$beta$beta$est),ncol(m$Par$beta$beta$est))
      else if(is.momentuHMM(m) & !is.null(m$mle$beta)) m$conditions$betaCons <- matrix(1:length(m$mle$beta),nrow(m$mle$beta),ncol(m$mle$beta))
    }
    if(is.null(m$conditions$betaRef)) m$conditions$betaRef <- as.integer(1:length(m$stateNames))
    if(is.momentuHMM(m)){
      if(is.null(m$mod$wpar)) m$mod$wpar <- m$mod$estimate
      if(is.null(m$mod$Sigma) & !is.null(m$mod$hessian)) m$mod$Sigma <- MASS::ginv(m$mod$hessian)
    } else {
      ####### compatability hack for change to MIcombine in momentuHMM >= 1.4.3 ######
      if(is.null(m$conditions$optInd)){
        for(i in names(m$conditions$dist)){
          m$conditions$cons[[i]]<-rep(1,length(m$conditions$cons[[i]]))
          m$conditions$workcons[[i]]<-rep(0,length(m$conditions$workcons[[i]]))
          m$conditions$workBounds[[i]]<-matrix(c(-Inf,Inf),nrow(m$conditions$workBounds[[i]]),2,byrow=TRUE)
        }
      }
      ################################################################################
    }
  } else if(!is.miHMM(m) & any(unlist(lapply(m,is.momentuHMM)))){
    m <- HMMfits(m)
  }
  m
}
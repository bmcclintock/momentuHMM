momentuHMMdists<-sort(c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm','norm','bern','vmConsensus','mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3','rw_norm','cat','negbinom','logis','t'))
moveHMMdists<-sort(c('gamma','weibull','exp','lnorm','wrpcauchy','vm'))
angledists<-sort(c('wrpcauchy','vm','vmConsensus'))
stepdists<-sort(c('gamma','weibull','exp','lnorm'))
singleParmdists<-sort(c('exp','pois','bern'))
nonnegativedists<-sort(c('gamma','weibull','exp','lnorm','pois','negbinom'))
zeroInflationdists<-sort(c('gamma','weibull','exp','lnorm','beta'))
oneInflationdists<-sort(c('beta'))
integerdists<-sort(c('bern','pois','cat','negbinom'))
mvndists <- c('mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3')
rwdists <- c('rw_norm','rw_mvnorm2','rw_mvnorm3')
splineList<-c("bs","ns","bSpline","mSpline","cSpline","iSpline")
meansList<-c("matrix","numeric","integer","logical","Date","POSIXlt","POSIXct","difftime")
meansListNoTime<-c("numeric","integer","logical")
plotArgs <- c("cex","cex.main","cex.lab","cex.axis","cex.legend","lwd","asp","legend.pos")
fitMethods<-c("nlm","Nelder-Mead","SANN")

#' @importFrom stats dbinom
dbern <- function (x, prob, log = FALSE) 
{
  return(stats::dbinom(x, 1, prob, log))
}

#' @importFrom stats pbinom
pbern <- function (q, prob, lower.tail = TRUE, log.p = FALSE) 
{
  return(stats::pbinom(q, 1, prob, lower.tail, log.p))
}

#' @importFrom stats rbinom
rbern <- function (n, prob) 
{
  return(stats::rbinom(n, 1, prob))
}

#' @importFrom stats dnbinom
dnegbinom <- function (x, mu, size, log = FALSE) 
{
  return(stats::dnbinom(x, size = size, mu = mu, log = log))
}

#' @importFrom stats pnbinom
pnegbinom <- function (q, mu, size, lower.tail = TRUE, log.p = FALSE) 
{
  return(stats::pnbinom(q, size = size, mu = mu, lower.tail = lower.tail, log.p = log.p))
}

#' @importFrom stats rnbinom
rnegbinom<- function (n, mu, size) 
{
  return(stats::rnbinom(n, size = size, mu = mu))
}

dmvnorm2 <- dmvnorm3 <- drw_mvnorm2 <- drw_mvnorm3 <- function(x,mean,sigma){
  dmvnorm_rcpp(x,mean,sigma)
}

rmvnorm2 <- rmvnorm3 <- rrw_mvnorm2 <- rrw_mvnorm3 <- mvtnorm::rmvnorm

drw_norm <- stats::dnorm
rrw_norm <- stats::rnorm

RWdata <- function(dist,data,knownStates){
  distnames <- names(dist)
  if(any(unlist(dist) %in% rwdists)){
    newdata <- NULL
    colInd <- NULL
    if(length(knownStates)){
      if("knownStates" %in% colnames(data)) stop("data cannot include a column named 'knownStates'")
      data$knownStates <- knownStates
    }
    ID <- unique(data$ID)
    for(j in ID){
      jInd <- which(data$ID==j)
      for(i in distnames){
        if(dist[[i]] %in% rwdists){
          tmpdata <- ldata <- data[jInd,,drop=FALSE]
          lInd <- 1:nrow(tmpdata)
          if(dist[[i]] %in% mvndists){
            if(inherits(data,"hierarchical")){
              iLevel <- attr(data,"coordLevel")
              lInd <- which(tmpdata$level==iLevel)
              ldata <- tmpdata[lInd,]
              colInd <- NULL
            }
            tmpdata[[paste0(i,".x_tm1")]] <- tmpdata[[paste0(i,".x")]]
            tmpdata[[paste0(i,".y_tm1")]] <- tmpdata[[paste0(i,".y")]]
            ldata[[paste0(i,".x_tm1")]] <- ldata[[paste0(i,".x")]]
            ldata[[paste0(i,".y_tm1")]] <- ldata[[paste0(i,".y")]]
            if(dist[[i]]=="rw_mvnorm2"){
              colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(distnames[which(!(distnames %in% i))],paste0(i,".x"),paste0(i,".y"))))]))
            } else if(dist[[i]]=="rw_mvnorm3"){
              tmpdata[[paste0(i,".z_tm1")]] <- tmpdata[[paste0(i,".z")]]
              ldata[[paste0(i,".z_tm1")]] <- ldata[[paste0(i,".z")]]
              colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(distnames[which(!(distnames %in% i))],paste0(i,".x"),paste0(i,".y"),paste0(i,".z"))))]))
            }
          } else {
            if(inherits(data,"hierarchical")){
              iLevel <- attr(data,"coordLevel")
              lInd <- which(tmpdata$level==iLevel)
              ldata <- tmpdata[lInd,]
              colInd <- NULL
            }
            tmpdata[[paste0(i,"_tm1")]] <- tmpdata[[i]]
            ldata[[paste0(i,"_tm1")]][lInd] <- ldata[[i]]
            colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% distnames))]))
          }
          if(inherits(data,"hierarchical")){
            ldata[,colInd] <- rbind(rep(NA,length(colInd)),ldata[-nrow(ldata),colInd])
            ldata <- ldata[-1,,drop=FALSE]
            tmpdata[lInd,colInd] <- rbind(rep(NA,length(colInd)),tmpdata[lInd[-length(lInd)],colInd])
            tmpdata <- tmpdata[-lInd[1],,drop=FALSE]
            tmpdata[lInd[-1]-1,colnames(tmpdata)] <- ldata[,colnames(tmpdata)]
            tmpdata[which(tmpdata$level!=iLevel),paste0(colnames(tmpdata)[!colnames(tmpdata) %in% colInd],"_tm1")] <- 0 # can't have NAs in covariates
          }
        }
      }
      if(!inherits(data,"hierarchical")){
        tmpdata[,colInd] <- rbind(rep(NA,length(colInd)),tmpdata[-nrow(tmpdata),colInd])
        tmpdata <- tmpdata[-1,,drop=FALSE]
      }
      newdata <- rbind(newdata,tmpdata)
    }
    class(newdata) <- class(data)
  } else newdata <- data
  newdata
}

#' @importFrom dplyr lag
crw <- function(x_tm1,lag=1){
  dplyr::lag(x_tm1,n=lag-1,default=x_tm1[1])-dplyr::lag(x_tm1,n=lag,default=x_tm1[1])
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
# this function maintains backwards compatibility with momentuHMM versions <1.4.0 (workBounds), <1.4.3 (betaCons), and <1.5.0 (mixtures)
delta_bc <- function(m){
  
  if(is.momentuHMM(m) | is.miSum(m)){
    if(!is.null(m$conditions$fit)){
      if(!m$conditions$fit) warning("The given model hasn't been fitted.")
    } else m$conditions$fit <- TRUE
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
    if(length(m$stateNames)>1 && is.null(m$conditions$betaCons)){
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
    if(is.null(m$conditions$mixtures)) m$conditions$mixtures <- 1
    if(is.null(m$covsPi)) m$covsPi <- matrix(1,length(unique(m$data$ID)),1)
    if(is.null(attr(m$data,"coords")) & !is.null(m$data$x) & !is.null(m$data$y)) attr(m$data,"coords") <- c("x","y")
  } else if(!is.miHMM(m) & any(unlist(lapply(m,is.momentuHMM)))){
    m <- HMMfits(m)
  }
  m
}
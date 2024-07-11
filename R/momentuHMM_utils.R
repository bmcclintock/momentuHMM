momentuHMMdists<-sort(c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm','norm','bern','vmConsensus','mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3','rw_norm','cat','negbinom','logis','t','ctds'))
moveHMMdists<-sort(c('gamma','weibull','exp','lnorm','wrpcauchy','vm'))
angledists<-sort(c('wrpcauchy','vm','vmConsensus'))
stepdists<-sort(c('gamma','weibull','exp','lnorm'))
singleParmdists<-sort(c('exp','pois','bern'))
nonnegativedists<-sort(c('gamma','weibull','exp','lnorm','pois','negbinom'))
zeroInflationdists<-sort(c('gamma','weibull','exp','lnorm','beta'))
oneInflationdists<-sort(c('beta'))
integerdists<-sort(c('bern','pois','cat','negbinom','ctds'))
mvndists <- c('mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3')
rwdists <- c('rw_norm','rw_mvnorm2','rw_mvnorm3')
CTHMMdists <- c(momentuHMMdists[which(!momentuHMMdists %in% c(angledists))])
CTHMMdepdists <- c(rwdists, "ctds","pois") # CTHMM distributions where observations and parameters depend on dt
splineList<-c("bs","ns",paste0("splines::",c("bs","ns")),"bSpline","mSpline","cSpline","iSpline",paste0("splines2::",c("bSpline","mSpline","cSpline","iSpline")),"poly","stats::poly")
meansList<-c("matrix","numeric","integer","Date","POSIXlt","POSIXct","difftime")
meansListNoTime<-c("numeric","integer")
modeList<- c("factor","logical")
plotArgs <- c("cex","cex.main","cex.lab","cex.axis","cex.legend","lwd","asp","legend.pos")
fitMethods<-c("nlm","TMB","Nelder-Mead","SANN")
badNames <- c("beta", "delta", "pi", "g0", "theta")

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

convertSigma <- function(sigcor,k){
  sigma <- sigcor
  for(i in 1:k){
    sigma[i,i] <- sigcor[i,i] * sigcor[i,i]
    for(j in 1:k){
      if(i!=j){
        sigma[i,j] <- sigcor[i,i] * sigcor[j,j] * sigcor[i,j]
        sigma[j,i] <- sigcor[i,i] * sigcor[j,j] * sigcor[i,j]
      }
    }
  }
  return(sigma)
}

rmvnorm2 <- rmvnorm3 <- rrw_mvnorm2 <- rrw_mvnorm3 <- function(n,mean,sigcor){
  mvtnorm::rmvnorm(n,mean,convertSigma(sigcor,length(mean)))
}

drw_norm <- stats::dnorm
rrw_norm <- stats::rnorm
prw_norm <- stats::pnorm

RWdata <- function(dist,data,knownStates,...){
  distnames <- names(dist)
  CT <- isTRUE(list(...)$CT)
  Time.name <- list(...)$Time.name
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
          if(inherits(data,"hierarchical")){
            iLevel <- attr(data,"coordLevel")
            lInd <- which(tmpdata$level==iLevel)
            ldata <- tmpdata[lInd,]
            colInd <- NULL
          }
          if(dist[[i]] %in% mvndists){
            tmpdata[[paste0(i,".x_tm1")]] <- tmpdata[[paste0(i,".x")]]
            tmpdata[[paste0(i,".y_tm1")]] <- tmpdata[[paste0(i,".y")]]
            ldata[[paste0(i,".x_tm1")]] <- ldata[[paste0(i,".x")]]
            ldata[[paste0(i,".y_tm1")]] <- ldata[[paste0(i,".y")]]
            if(dist[[i]]=="rw_mvnorm2"){
              #if(!CT) 
              colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(paste0(i,".x"),paste0(i,".y"))))]))
              ##else if(inherits(data,"hierarchical")) colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(Time.name,"dt",paste0(i,".x"),paste0(i,".y"))))]))
              #else colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(Time.name,paste0(i,".x"),paste0(i,".y"))))]))
            } else if(dist[[i]]=="rw_mvnorm3"){
              tmpdata[[paste0(i,".z_tm1")]] <- tmpdata[[paste0(i,".z")]]
              ldata[[paste0(i,".z_tm1")]] <- ldata[[paste0(i,".z")]]
              #if(!CT) 
              colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(paste0(i,".x"),paste0(i,".y"),paste0(i,".z"))))]))
              ##else if(inherits(data,"hierarchical")) colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(Time.name,"dt",paste0(i,".x"),paste0(i,".y"),paste0(i,".z"))))]))
              #else colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(Time.name,paste0(i,".x"),paste0(i,".y"),paste0(i,".z"))))]))
            }
          } else {
            tmpdata[[paste0(i,"_tm1")]] <- tmpdata[[i]]
            ldata[[paste0(i,"_tm1")]][lInd] <- ldata[[i]]
            #if(!CT) 
            colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% i))]))
            ##else if(inherits(data,"hierarchical")) colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(distnames,Time.name,"dt")))]))
            #else colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(distnames,Time.name)))]))
          }
          if(inherits(data,"hierarchical")){
            ##if(CT) tmpdata$dt[lInd[1]] <- 0
            ldata[,colInd] <- rbind(rep(NA,length(colInd)),ldata[-nrow(ldata),colInd])
            ldata <- ldata[-1,,drop=FALSE]
            tmpdata[lInd[-length(lInd)],colnames(tmpdata)] <- ldata[,colnames(tmpdata)]
            tmpdata[lInd[length(lInd)],colnames(tmpdata)[which(!(colnames(tmpdata) %in% colInd))]] <- NA
            tmpdata[which(tmpdata$level!=iLevel),paste0(colnames(tmpdata)[!colnames(tmpdata) %in% c(Time.name,"dt",colInd,distnames[which(!(distnames %in% i))])],"_tm1")] <- 0 # can't have NAs in covariates
          }
        }
      }
      if(!inherits(data,"hierarchical")){
        lastRow <- tmpdata[nrow(tmpdata),]
        tmpdata[,colInd] <- rbind(rep(NA,length(colInd)),tmpdata[-nrow(tmpdata),colInd])
        tmpdata <- tmpdata[-1,,drop=FALSE]
        tmpdata <- rbind(tmpdata,lastRow)
        tmpdata[nrow(tmpdata),colnames(tmpdata)[which(!(colnames(tmpdata) %in% colInd))]] <- NA
      }
      newdata <- rbind(newdata,tmpdata)
    }
    class(newdata) <- class(data)
  } else newdata <- data
  newdata
}

# @importFrom dplyr lag
crw <- function(x_tm1,lag=1,dt){
  for(pkg in c("dplyr")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package \"",pkg,"\" needed for crw function to work. Please install it.",
           call. = FALSE)
    }
  }
  if(missing(dt)) dt <- rep(1,length(x_tm1))
  else if(length(dt)==1) dt <- rep(dt,length(x_tm1))
  else if(length(dt)!=length(x_tm1)) stop("argument 'dt' in 'crw' function must be of length 1 or ",length(x_tm1))
  
  diff <- dplyr::lag(x_tm1,n=lag-1,default=x_tm1[1])-dplyr::lag(x_tm1,n=lag,default=x_tm1[1])
  vel <- diff/dplyr::lag(dt,n=lag,default=0)
  if(any(!is.finite(vel))){
    vel[which(!is.finite(vel))] <- 0
  }
  return(vel)
}

radian <- function(degree) 
{
  radian <- degree * (pi/180)
  radian
}

ctPar <- function(par,dist,nbStates,data){
  for(i in names(dist)){
    if(dist[[i]] %in% rwdists){
      par[[i]][1:nbStates,] <- t(apply(par[[i]][1:nbStates,,drop=FALSE] - rep(data[[paste0(i,".x_tm1")]],each=nbStates),1,function(x) x*data$dt)) + rep(data[[paste0(i,".x_tm1")]],each=nbStates)
      par[[i]][nbStates+1:nbStates,] <- t(apply(par[[i]][nbStates+1:nbStates,,drop=FALSE] - rep(data[[paste0(i,".y_tm1")]],each=nbStates),1,function(x) x*data$dt)) + rep(data[[paste0(i,".y_tm1")]],each=nbStates)
      if(dist[[i]]=="rw_mvnorm3"){
        par[[i]][2*nbStates+1:nbStates,] <- t(apply(par[[i]][2*nbStates+1:nbStates,,drop=FALSE] - rep(data[[paste0(i,".z_tm1")]],each=nbStates),1,function(x) x*data$dt)) + rep(data[[paste0(i,".z_tm1")]],each=nbStates)
        par[[i]][3*nbStates+1:(3*nbStates),] <- t(apply(par[[i]][3*nbStates+1:(3*nbStates),,drop=FALSE],1,function(x) x*sqrt(data$dt)))
      } else {
        par[[i]][2*nbStates+1:(2*nbStates),] <- t(apply(par[[i]][2*nbStates+1:(2*nbStates),,drop=FALSE],1,function(x) x*sqrt(data$dt)))
      }
    } else if(dist[[i]]=="ctds"){
      p <- matrix(0,nrow(par[[i]])+nbStates,ncol(par[[i]]))
      par[[i]] <- matrix(t(apply(log(par[[i]]),1,function(x) x+log(data$dt))),ncol=ncol(par[[i]]))
      for(j in 1:nbStates){
        catInd <- seq(j,nrow(par[[i]]),nbStates)
        probPar <- rbind(par[[i]][catInd,,drop=FALSE],rep(0,ncol(par[[i]])))
        expPar <- exp(probPar)
        prob <- expPar/rep(colSums(expPar),each=length(catInd)+1)
        for(k in which(!is.finite(colSums(prob)))){
          tmp <- exp(Brobdingnag::as.brob(probPar[,k]))
          prob[,k] <- as.numeric(tmp/Brobdingnag::sum(tmp))
        }
        p[seq(j,nrow(par[[i]])+nbStates,nbStates),] <- prob
      }
      par[[i]] <- p
    } else if(dist[[i]] %in% CTHMMdepdists) {
      par[[i]] <- matrix(t(apply(par[[i]],1,function(x) x*data$dt)),ncol=ncol(par[[i]]))
      #par[[i]] <- t(apply(par[[i]],1,function(x) x*data$dt))
      #par[[i]] <- matrix(apply(par[[i]],1,function(x) x*data$dt),ncol=ncol(par[[i]]))
    }
  }
  par
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
  if(any(grepl("Foreach loop \\(doParallelMC\\) had changed the current RNG type: RNG was restored to same type, next state",w))
     | any(grepl("already exporting variable\\(s\\):",w)))
    invokeRestart("muffleWarning")
  else w
}

# suppress CT argument warning
muffleCTwarning <- function(w) {
  q <- iconv(c(intToUtf8(8216), intToUtf8(8217)),"UTF-8", "")
  if(any(grepl(paste0("extra argument ",q[1],"CT",q[2]," will be disregarded"),w)) |
     any(grepl(paste0("extra arguments ",q[1],"CT",q[2],", ",
                      q[1],"kappa",q[2]," will be disregarded"),w)) |
     any(grepl(paste0("extra arguments ",q[1],"CT",q[2],", ",
                      q[1],"Time.name",q[2],", ",
                      q[1],"kappa",q[2]," will be disregarded"),w)) |
     any(grepl(paste0("extra arguments ",q[1],"CT",q[2],", ",
                      q[1],"Time.name",q[2],", ",
                      q[1],"kappa",q[2]," will be disregarded"),w)) |
     any(grepl(paste0("extra arguments ",q[1],"CT",q[2],", ",
                      q[1],"Time.name",q[2],", ",
                      q[1],"kappa",q[2],", ",
                      q[1],"dontcheck",q[2]," will be disregarded"),w)) |
     any(grepl(paste0("extra argument ",q[1],"dontcheck",q[2]," will be disregarded"),w)))
    invokeRestart("muffleWarning")
  else w
}

muffleCTDSwarning <- function(w) {
  q <- iconv(c(intToUtf8(8216), intToUtf8(8217)),"UTF-8", "")
  if(any(grepl(paste0("extra arguments ",q[1],"CT",q[2],", ",
                                         q[1],"kappa",q[2],", ",
                                         q[1],"ctds",q[2],", ",
                                         q[1],"rast",q[2],", ",
                                         q[1],"directions",q[2],", ",
                                         q[1],"moveState",q[2]," will be disregarded"),w)))
    invokeRestart("muffleWarning")
  else w
}

chkDotsCT <- function(...){
  if("CT" %in% names(list(...)) && !any(grepl("CTHMM",unlist(lapply(sys.calls(),function(x) deparse(x)[1]))) | 
                                        grepl("CTDS",unlist(lapply(sys.calls(),function(x) deparse(x)[1]))))) stop(sprintf("In %s :\n extra argument 'CT' is invalid", 
                                                                                                                              paste(deparse(sys.call()[[1]], control = c()), 
                                                                                                                                    collapse = "\n")), call. = FALSE, domain = NA)
  if("ctds" %in% names(list(...)) && !any(grepl("CTDS",unlist(lapply(sys.calls(),function(x) deparse(x)[1]))))) stop(sprintf("In %s :\n extra argument 'ctds' is invalid", 
                                                                                                                                paste(deparse(sys.call()[[1]], control = c()), 
                                                                                                                                      collapse = "\n")), call. = FALSE, domain = NA)
  if("rast" %in% names(list(...)) && !any(grepl("CTDS",unlist(lapply(sys.calls(),function(x) deparse(x)[1]))))) stop(sprintf("In %s :\n extra argument 'rast' is invalid", 
                                                                                                                                paste(deparse(sys.call()[[1]], control = c()), 
                                                                                                                                      collapse = "\n")), call. = FALSE, domain = NA) 
  if("directions" %in% names(list(...)) && !any(grepl("CTDS",unlist(lapply(sys.calls(),function(x) deparse(x)[1]))))) stop(sprintf("In %s :\n extra argument 'directions' is invalid", 
                                                                                                                                      paste(deparse(sys.call()[[1]], control = c()), 
                                                                                                                                            collapse = "\n")), call. = FALSE, domain = NA) 
}

noMove <- function(out,directions){
  nomove <- factor(,levels=c(FALSE,TRUE))
  for(i in unique(out$ID)){
    iInd <- which(out$ID==i)
    nomove[iInd[1]] <- FALSE
    nomove[iInd[-1]] <- factor(out$z[iInd[1:(length(iInd)-1)]]==(directions+1))
    
  }
  return(nomove)
}

# .combine function for multiple rbinds in foreach
comb <- function(x, ...) {  
    mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# #' @importFrom doFuture registerDoFuture
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar% foreach
# #' @importFrom future multisession plan
# #' @importFrom iterators icount
progBar <- function(kk, N, per = 1) {
  if (kk %in% seq(0, N, per)) {
    x <- round(kk * 100 / N)
    message("[ ", 
            paste(rep("=", x), collapse = ""),
            paste(rep("-", 100 - x), collapse = ""), 
            " ] ", x, "%", "\r",
            appendLF = FALSE)
    if (kk == N) message("\n")
  }
}

installDataTree <- function(){
  if (!requireNamespace("data.tree", quietly = TRUE)) {
    stop("Package \"data.tree\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
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
          m$conditions$workBounds[[i]]<-matrix(c(-Inf,Inf),nrow(m$conditions$workBounds[[i]]),2,byrow=TRUE)
        }
      }
      ################################################################################
    }
    if(is.null(m$conditions$mixtures)) m$conditions$mixtures <- 1
    if(is.null(m$covsPi)) m$covsPi <- matrix(1,length(unique(m$data$ID)),1)
    if(is.null(attr(m$data,"coords")) & !is.null(m$data$x) & !is.null(m$data$y)) attr(m$data,"coords") <- c("x","y")
    if(is.null(attr(m$data,"gradient"))) attr(m$data,"gradient") <- FALSE
    if(is.null(m$conditions$kappa)) m$conditions$kappa <- Inf
  } else if(!is.miHMM(m) & any(unlist(lapply(m,is.momentuHMM)))){
    m <- HMMfits(m)
  }
  m
}

#' Transform a raster into a (x,y,z) list
#'
#' @param rast \code{\link[raster]{raster}} layer object for spatially referenced covariates.
#' @return List of three elements: x (grid of x values), y (grid of y values), and z (matrix of values).
#' @export
collapseRaster <- function(rast){
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package \"raster\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  lim <- as.vector(raster::extent(rast))
  rres <- raster::res(rast)
  xgrid <- seq(lim[1] + rres[1]/2, lim[2] - rres[1]/2, by = rres[1])
  ygrid <- seq(lim[3] + rres[2]/2, lim[4] - rres[2]/2, by = rres[2])
  z <- t(apply(raster::as.matrix(rast), 2, rev))
  return(list(x = xgrid, y = ygrid, z = z))
}

gridCell <- function(loc, xgrid, ygrid, covmat){
  ix <- findInterval(loc[1], xgrid)
  iy <- findInterval(loc[2], ygrid)
  coords <- c(xgrid[ix], xgrid[ix + 1], ygrid[iy], ygrid[iy + 1])
  values <- covmat[ix:(ix + 1), iy:(iy + 1)]
  return(list(coords = coords, values = values))
}

biGrad <- function(loc, covlist) {
  J <- length(covlist)
  grad_val <- sapply(1:J, function(j) {
    x_grid <- covlist[[j]]$x
    y_grid <- covlist[[j]]$y
    covmat <- covlist[[j]]$z
    cell <- tryCatch(gridCell(loc = loc, xgrid = x_grid, ygrid = y_grid, 
                     covmat = covmat),error=function(e) e)
    if(!inherits(cell,"error") && length(cell$coords)==4 && all(dim(cell$values)==2)){
      x <- cell$coords[1:2]
      y <- cell$coords[3:4]
      f <- cell$values
      dfdx <- ((y[2] - loc[2]) * (f[2, 1] - f[1, 1]) + (loc[2] - 
                                                          y[1]) * (f[2, 2] - f[1, 2]))/((y[2] - y[1]) * (x[2] - 
                                                                                                           x[1]))
      dfdy <- ((x[2] - loc[1]) * (f[1, 2] - f[1, 1]) + (loc[1] - 
                                                          x[1]) * (f[2, 2] - f[2, 1]))/((y[2] - y[1]) * (x[2] - 
                                                                                                           x[1]))
      return(c(dfdx, dfdy))
    } else return(c(0,0))
  })
  return(grad_val)
}

biGradArray <- function(locs, covlist){
  if (!inherits(locs, "matrix")) 
    stop("'locs' must be a matrix")
  grad <- unlist(lapply(1:nrow(locs), function(i) biGrad(loc = locs[i,], covlist = covlist)))
  gradarray <- array(grad, c(2, length(covlist), nrow(locs)))
  gradarray <- aperm(gradarray, c(3, 1, 2))
  return(gradarray)
}

#' Calculate gradient of spatial covariates using bilinear interpolation
#'
#' @param data Data frame of data streams. At a minimum, it must contain fields matching \code{coordNames}
#' @param spatialCovs List of \code{\link[raster]{raster}} layer objects for spatially referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data for each time step.
#' @param collapseRast List of collapsed \code{\link[raster]{raster}} layer objects (see \code{\link{collapseRaster}}). Ignored unless \code{spatialCovs} is missing.
#' @param coordNames Names of the coordinates in \code{data}. Default: \code{c("x","y")}.
#' @return The gradients are appended to \code{data} with ``\code{.x}'' (easting gradient) and ``\code{.y}'' (northing gradient) suffixes added to the names of \code{spatialCovs}. For example, if \code{cov1} is the name of a spatial covariate, then the returned \code{data} object will include the fields ``\code{cov1.x}'' and ``\code{cov1.y}''.
#' @export
getGradients <- function(data, spatialCovs, collapseRast, coordNames=c("x","y")){
  if(!missing(spatialCovs)) covlist <- lapply(spatialCovs, collapseRaster)
  else covlist <- collapseRast
  gradarray <- biGradArray(locs = as.matrix(data[coordNames]), covlist = covlist)
  covNames <- names(covlist)
  gradNames <- c("x","y")
  for(i in 1:length(covNames)){
    for(j in 1:length(gradNames)){
      data[[paste0(covNames[i],".",gradNames[j])]] <- gradarray[,j,i]
    }
  }
  return(data)
}

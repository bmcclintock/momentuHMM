# format miSum object for use in momentuHMM object functions
formatmiSum <- function(miSum){
  miSum$mle <- lapply(miSum$Par$real,function(x) x$est)
  miSum$mle$beta <- miSum$Par$beta$beta$est
  miSum$mle$pi <- miSum$Par$real$pi$est
  miSum$mle$delta <- miSum$Par$real$delta$est
  miSum$mod <- list()
  if(!is.null(miSum$conditions$recharge)){
    nbRecovs <- ncol(model.matrix(miSum$conditions$recharge$g0,miSum$g0covs))-1
    nbG0covs <- ncol(model.matrix(miSum$conditions$recharge$theta,miSum$reCovs))-1
    miSum$mle$g0 <- c(miSum$Par$beta$g0$est)
    names(miSum$mle$g0) <- colnames(miSum$Par$beta$g0$est)
    miSum$mle$theta <- c(miSum$Par$beta$theta$est)
    names(miSum$mle$theta) <- colnames(miSum$Par$beta$theta$est)
  } else nbRecovs <- nbG0covs <- 0
  miSum$mod$estimate <- expandPar(miSum$MIcombine$coefficients,miSum$conditions$optInd,unlist(miSum$conditions$fixPar),miSum$conditions$wparIndex,miSum$conditions$betaCons,miSum$conditions$deltaCons,length(miSum$stateNames),ncol(miSum$covsDelta)-1,miSum$conditions$stationary,nrow(miSum$Par$beta$beta$est)/miSum$conditions$mixtures-1,nbRecovs+nbG0covs,miSum$conditions$mixtures,ncol(miSum$covsPi)-1)
  miSum$mod$wpar <- miSum$MIcombine$coefficients
  miSum
}
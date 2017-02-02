#' @export
#' @importFrom crawl crwPostIS crwSimulator
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
MIfitHMM<-function(nSims,ncores,crwOut,
                   nbStates, dist, Par, beta0 = NULL, delta0 = NULL,
                   estAngleMean = NULL, formula = ~1, stationary = FALSE, verbose = 0,
                   nlmPar = NULL, fit = TRUE, DM = NULL, cons = NULL,
                   userBounds = NULL, workcons = NULL, stateNames = NULL, knownStates=NULL,
                   type=c("LL", "UTM"),coordNames=c("x","y"),covNames=NULL,spatialCovs=NULL,
                   method = "IS", parIS = 1000, dfSim = Inf, grid.eps = 1, crit = 2.5, scaleSim = 1, force.quad,
                   fullPost = TRUE, dfPostIS = Inf, scalePostIS = 1,thetaSamp = NULL,
                   poolEstimates = TRUE){
  
  if(!is.crwData(crwOut))
    stop("crwOut must be a list containing 'crwFits' (a list of crwFit objects) and 'crwPredict' (a crwPredict object)")
  
  Time.name<-attr(crwOut$crwPredict,"Time.name")
  distnames<-names(dist)[which(!(names(dist) %in% c("step","angle")))]
  ids = unique(crwOut$crwPredict$ID)
  
  if(nSims>1){
    cat('Drawing',nSims,'realizations from the position process using crawl::crwPostIs...')
    
    registerDoParallel(cores=ncores)
    crwSim <- foreach(i = 1:length(ids), .export="crwSimulator") %dopar% {
      crawl::crwSimulator(crwOut$crwFits[[i]],predTime=crwOut$crwPredict[[Time.name]][which(crwOut$crwPredict$ID==ids[i] & crwOut$crwPredict$locType=="p")], method = method, parIS = parIS,
                                df = dfSim, grid.eps = grid.eps, crit = crit, scale = scaleSim, force.quad)
    }
    stopImplicitCluster()
    
    registerDoParallel(cores=ncores)
    mh<-
      foreach(j = 1:nSims, .export=c("crwPostIS","prepData")) %dopar% {
        locs<-data.frame()
        for(i in 1:length(ids)){
          tmp<-tryCatch({crawl::crwPostIS(crwSim[[i]], fullPost = fullPost, df = dfPostIS, scale = scalePostIS, thetaSamp = thetaSamp)},error=function(e) e)
          if(!all(class(tmp) %in% c("crwIS","list"))) stop('crawl::crwPostIS error for individual ',ids[i],'; ',tmp,'  Check crwPostIS arguments, crawl::crwMLE model fits, and/or consult crawl documentation.')
          locs<-rbind(locs,tmp$alpha.sim[,c("mu.x","mu.y")])
        }
        df<-data.frame(x=locs$mu.x,y=locs$mu.y,crwOut$crwPredict[,c("ID",distnames,covNames),drop=FALSE])[which(crwOut$crwPredict$locType=="p"),]
        prepData(df,type=type,coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs)
      }
    stopImplicitCluster()
    cat("DONE\n")
    cat('Fitting',nSims,'realizations of the position process using fitHMM...')
  } else {
    mh <- list()
    df <- data.frame(x=crwOut$crwPredict$mu.x,y=crwOut$crwPredict$mu.y,crwOut$crwPredict[,c("ID",distnames,covNames),drop=FALSE])[which(crwOut$crwPredict$locType=="p"),]
    mh[[1]] <- prepData(df,type=type,coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs)
    cat('Fitting the most likely position process using fitHMM...')
  }

  registerDoParallel(cores=ncores)
  fits <-
    foreach(j = 1:nSims, .export=c("fitHMM"), .errorhandling="pass") %dopar% {

      fitHMM(mh[[j]],nbStates, dist, Par, beta0, delta0,
                              estAngleMean, formula, stationary, verbose,
                              nlmPar, fit, DM, cons,
                              userBounds, workcons, stateNames, knownStates)
    }  
  stopImplicitCluster()
  cat("DONE\n")
  
  for(i in which(!unlist(lapply(fits,function(x) inherits(x,"momentuHMM"))))){
    warning('Fit #',i,' failed; ',fits[[i]])
  }
  
  if(nSims==1) out<-fits[[1]]
  else {
    if(poolEstimates) out<-MI_summary(fits)
    else out <- fits
  }
  out
}
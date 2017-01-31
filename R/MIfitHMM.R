#' @export
#' @importFrom crawl crwPostIS
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
MIfitHMM<-function(nSims,ncores,obsData,crwFits,
                   nbStates, dist, Par, beta0 = NULL, delta0 = NULL,
                   estAngleMean = NULL, formula = ~1, stationary = FALSE, verbose = 0,
                   nlmPar = NULL, fit = TRUE, DM = NULL, cons = NULL,
                   userBounds = NULL, workcons = NULL, stateNames = NULL, knownStates=NULL,
                   type=c("LL", "UTM"),coordNames=c("x","y"),covNames=NULL,spatialCovs=NULL,
                   fullPost = TRUE, df = Inf, scale = 1,thetaSamp = NULL){
  
  if(nrow(obsData)!=nrow(crwFits$predData)){
    if(!length(names(dist)[which(!(names(dist) %in% c("step","angle")))])){
      if(nrow(obsData)==sum(!is.na(crwFits$predData$x))){
        Time.name<-crwFits$crwSim[[1]]$Time.name
        tmpObs <- data.frame(ID=crwFits$predData$ID,crwFits$predData[Time.name])
        obsData<-merge(obsData,tmpObs,by=c("ID",Time.name),all.y=TRUE)
      } else stop('obsData should consist of ',sum(!is.na(crwFits$predData$x)),' rows')
    } else {
        stop('Because additional data streams (',paste0(names(dist)[which(!(names(dist) %in% c("step","angle")))],collapse=", "),') are used, obsData must include ',nrow(crwFits$predData),' rows. 
        These should include ',sum(crwFits$predData$locType=="o"),' rows corresponding to each observed location, 
        and ',sum(crwFits$predData$locType=="p"),' rows corresponding to the predicted time steps (sorted by ID and time).  See crwFits$predData for appropriate structure.')
    }
  }
  
  nbAnimals<-length(crwFits$model_fits)
  distnames<-names(dist)[which(!(names(dist) %in% c("step","angle")))]
  
  if(nSims>1){
    cat('Drawing',nSims,'realizations from the position process using crawl::crwPostIs...')
    registerDoParallel(cores=ncores)
    mh<-
      foreach(j = 1:nSims, .export=c("crwPostIS","prepData")) %dopar% {
        locs<-data.frame()
        locType<-character()
        for(i in 1:nbAnimals){
          tmp<-tryCatch({crawl::crwPostIS(crwFits$crwSim[[i]], fullPost, df = df, scale = scale, thetaSamp = thetaSamp)},error=function(e) e)
          if(!all(class(tmp) %in% c("crwIS","list"))) stop('crawl::crwPostIS error for individual ',i,'; ',tmp,'  Check crwPostIS arguments, crwFits$model_fits, and/or consult crawl documentation.')
          locs<-rbind(locs,tmp$alpha.sim[,c("mu.x","mu.y")])
          locType<-c(locType,tmp$locType)
        }
        predData<-list(mu.x=locs$mu.x,mu.y=locs$mu.y,locType=locType)
        df<-data.frame(x=predData$mu.x,y=predData$mu.y,obsData[,c("ID",distnames),drop=FALSE])[which(predData$locType=="p"),]
        prepData(df,type=type,coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs)
      }
    stopImplicitCluster()
    cat("DONE\n")
    cat('Fitting',nSims,'realizations of the position process using fitHMM...')
  } else {
    mh <- list()
    df <- data.frame(x=crwFits$predData$mu.x,y=crwFits$predData$mu.y,obsData[,c("ID",distnames),drop=FALSE])[which(crwFits$predData$locType=="p"),]
    mh[[1]] <- prepData(df,type=type,coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs)
    cat('Fitting the most likely position process using fitHMM...')
  }

  registerDoParallel(cores=ncores)
  fits <-
    foreach(j = 1:nSims, .export=c("fitHMM")) %dopar% {

      fit<-fitHMM(mh[[j]],nbStates, dist, Par, beta0, delta0,
                              estAngleMean, formula, stationary, verbose,
                              nlmPar, fit, DM, cons,
                              userBounds, workcons, stateNames, knownStates)
    }  
  stopImplicitCluster()
  cat("DONE\n")
  
  if(nSims==1) fits<-fits[[1]]
  
  fits
}
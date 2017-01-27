#' @export
#' @importFrom crawl crwPostIS
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
MIfitHMM<-function(n,ncores,obsData,crwFits,
                   nbStates, dist, Par, beta0 = NULL, delta0 = NULL,
                   estAngleMean = NULL, formula = ~1, stationary = FALSE, verbose = 0,
                   nlmPar = NULL, fit = TRUE, DM = NULL, cons = NULL,
                   userBounds = NULL, workcons = NULL, stateNames = NULL,
                   type=c("LL", "UTM"),coordNames=c("x","y"),covNames=NULL,spatialCovs=NULL,
                   fullPost = TRUE, df = Inf, scale = 1,thetaSamp = NULL){
  
  nbAnimals<-length(crwFits$model_fits)
  
  registerDoParallel(cores=ncores)
  multipleFits <-
    foreach(j = 1:n, .export=c("crwPostIS","prepData","fitHMM")) %dopar% {
      locs<-data.frame()
      locType<-character()
      for(i in 1:nbAnimals){
        tmp<-crawl::crwPostIS(crwFits$crwSim[[i]], fullPost, df = df, scale = scale, thetaSamp = thetaSamp)
        locs<-rbind(locs,tmp$alpha.sim[,c("mu.x","mu.y")])
        locType<-c(locType,tmp$locType)
      }
      predData<-list(mu.x=locs$mu.x,mu.y=locs$mu.y,locType=locType)
      distnames<-names(dist)[which(!(names(dist) %in% c("step","angle")))]
      mh<-data.frame(x=predData$mu.x,y=predData$mu.y,obsData[,c("ID",distnames),drop=FALSE])[which(predData$locType=="p"),]
      mh<-prepData(mh,type=type,coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs)
      fit<-fitHMM(mh,nbStates, dist, Par, beta0, delta0,
                              estAngleMean, formula, stationary, verbose,
                              nlmPar, fit, DM, cons,
                              userBounds, workcons, stateNames)
    }  
  stopImplicitCluster()
  multipleFits
}
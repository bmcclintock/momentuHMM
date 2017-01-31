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
  
  registerDoParallel(cores=ncores)
  multipleFits <-
    foreach(j = 1:n, .export=c("crwPostIS","prepData","fitHMM")) %dopar% {
      predData<-tryCatch({
        locs<-data.frame()
        locType<-character()
        for(i in 1:nbAnimals){
          tmp<-crawl::crwPostIS(crwFits$crwSim[[i]], fullPost, df = df, scale = scale, thetaSamp = thetaSamp)
          locs<-rbind(locs,tmp$alpha.sim[,c("mu.x","mu.y")])
          locType<-c(locType,tmp$locType)
        }
        list(mu.x=locs$mu.x,mu.y=locs$mu.y,locType=locType)},error=function(e) "badSim")
      if(predData!="badSim"){
        distnames<-names(dist)[which(!(names(dist) %in% c("step","angle")))]
        mh<-data.frame(x=predData$mu.x,y=predData$mu.y,obsData[,c("ID",distnames),drop=FALSE])[which(predData$locType=="p"),]
        mh<-prepData(mh,type=type,coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs)
        fit<-fitHMM(mh,nbStates, dist, Par, beta0, delta0,
                                estAngleMean, formula, stationary, verbose,
                                nlmPar, fit, DM, cons,
                                userBounds, workcons, stateNames)
      } else fit<-predData
      fit
    }  
  stopImplicitCluster()
  multipleFits
}
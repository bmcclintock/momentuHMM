#' @export
#' @importFrom crawl crwPostIS crwSimulator
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
MIfitHMM<-function(nSims, ncores, poolEstimates = TRUE, alpha = 0.95,
                   miData,
                   nbStates, dist, Par, beta0 = NULL, delta0 = NULL,
                   estAngleMean = NULL, 
                   circularAngleMean = NULL,
                   formula = ~1, stationary = FALSE, verbose = 0,
                   nlmPar = NULL, fit = TRUE, DM = NULL, cons = NULL,
                   userBounds = NULL, workcons = NULL, stateNames = NULL, knownStates=NULL,
                   covNames=NULL,spatialCovs=NULL,centers=NULL,
                   method = "IS", parIS = 1000, dfSim = Inf, grid.eps = 1, crit = 2.5, scaleSim = 1, force.quad,
                   fullPost = TRUE, dfPostIS = Inf, scalePostIS = 1,thetaSamp = NULL
                   ){

  j <- NULL #gets rid of no visible binding for global variable ‘j’ NOTE in R cmd check
  
  if(is.crwData(miData)){
    
    model_fits <- miData$crwFits
    predData <- miData$crwPredict

    Time.name<-attr(predData,"Time.name")
    ids = unique(predData$ID)
    distnames<-names(dist)[which(!(names(dist) %in% c("step","angle")))]
    coordNames <- attr(predData,"coord")
    
    if(nSims>1){
      if(all(unlist(lapply(model_fits,function(x) is.null(x$err.model))))) stop("Multiple realizations of the position process cannot be drawn if there is no location measurement error")
      cat('Drawing',nSims,'realizations from the position process using crawl::crwPostIs...')
      
      registerDoParallel(cores=ncores)
      crwSim <- foreach(i = 1:length(ids), .export="crwSimulator") %dopar% {
        if(!is.null(model_fits[[i]]$err.model))
          crawl::crwSimulator(model_fits[[i]],predTime=predData[[Time.name]][which(predData$ID==ids[i] & predData$locType=="p")], method = method, parIS = parIS,
                                  df = dfSim, grid.eps = grid.eps, crit = crit, scale = scaleSim, force.quad)
      }
      stopImplicitCluster()
      
      registerDoParallel(cores=ncores)
      miData<-
        foreach(j = 1:nSims, .export=c("crwPostIS","prepData")) %dopar% {
          locs<-data.frame()
          for(i in 1:length(ids)){
            if(!is.null(model_fits[[i]]$err.model)){
              tmp<-tryCatch({crawl::crwPostIS(crwSim[[i]], fullPost = fullPost, df = dfPostIS, scale = scalePostIS, thetaSamp = thetaSamp)},error=function(e) e)
              if(!all(class(tmp) %in% c("crwIS","list"))) stop('crawl::crwPostIS error for individual ',ids[i],'; ',tmp,'  Check crwPostIS arguments, crawl::crwMLE model fits, and/or consult crawl documentation.')
              locs<-rbind(locs,tmp$alpha.sim[,c("mu.x","mu.y")])
            } else {
              locs<-rbind(locs,predData[which(predData$ID==ids[i]),c("mu.x","mu.y")])
            }
          }
          df<-data.frame(x=locs$mu.x,y=locs$mu.y,predData[,c("ID",distnames,covNames),drop=FALSE])[which(predData$locType=="p"),]
          prepData(df,type="UTM",coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs,centers=centers)
        }
      stopImplicitCluster()
      cat("DONE\n")
      cat('Fitting',nSims,'realizations of the position process using fitHMM...')
    } else {
      miData <- list()
      df <- data.frame(x=predData$mu.x,y=predData$mu.y,predData[,c("ID",distnames,covNames),drop=FALSE])[which(predData$locType=="p"),]
      miData[[1]] <- prepData(df,type="UTM",coordNames=coordNames,covNames=covNames,spatialCovs=spatialCovs,centers=centers)
      cat('Fitting the most likely position process using fitHMM...')
    }
    
  } else {
    if(!is.list(miData)) stop("miData must either be a crwData object (as returned by CRAWLwrap) or a list of momentuHMMData objects (as returned by simData and prepData)")
    if(!any(unlist(lapply(miData,function(x) inherits(x,"momentuHMMData"))))) stop("miData must either be a crwData object (as returned by CRAWLwrap) or a list of momentuHMMData objects (as returned by simData and prepData)")
    if(nSims>length(miData))
      stop("nSims is greater than the length of miData. nSims must be <=",length(miData))
    cat('Fitting',nSims,'imputation(s) using fitHMM...')
  }
  
  registerDoParallel(cores=ncores)
  fits <-
    foreach(j = 1:nSims, .export=c("fitHMM"), .errorhandling="pass") %dopar% {

      suppressMessages(fitHMM(miData[[j]],nbStates, dist, Par, beta0, delta0,
                              estAngleMean, circularAngleMean, formula, stationary, verbose,
                              nlmPar, fit, DM, cons,
                              userBounds, workcons, stateNames, knownStates))
    }  
  stopImplicitCluster()
  cat("DONE\n")
  
  for(i in which(!unlist(lapply(fits,function(x) inherits(x,"momentuHMM"))))){
    warning('Fit #',i,' failed; ',fits[[i]])
  }
  
  if(nSims==1) out<-fits[[1]]
  else {
    if(poolEstimates) out<-MI_summary(fits,alpha=alpha,ncores=ncores,includeHMMfits=TRUE)
    else out <- fits
  }
  out
}
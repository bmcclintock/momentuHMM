#' @export
#' @importFrom sp coordinates
#' @importFrom crawl crwMLE crwPredict crwSimulator
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats formula

CRAWLwrap<-function(obsData, timeStep=1, ncores, 
                    mov.model = NULL, err.model = NULL, activity = NULL, drift = NULL, 
                    coord = c("x", "y"), Time.name = "time", initial.state, theta, fixPar, 
                    methodMLE = "L-BFGS-B", control = NULL, constr = NULL, 
                    prior = NULL, need.hess = TRUE, initialSANN = list(maxit = 200), attempts = 1,
                    predTime = NULL, 
                    methodSim = "IS", parIS = 1000, df = Inf, grid.eps = 1, crit = 2.5, scale = 1, force.quad)
{
  
  if(any(!(c("ID",Time.name,coord) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name,coord)[!(c("ID",Time.name,coord) %in% names(obsData))],collapse=","))
  
  toProj = obsData[!is.na(obsData$y),c("ID",Time.name,coord)]
  sp::coordinates(toProj) = formula(paste0("~",paste0(coord,collapse="+")))
  
  ids = unique(obsData$ID)
  ind_data<-vector('list',length(ids))
  
  for(i in 1:length(ids)){
    ind_data[[i]] = obsData[which(obsData$ID==ids[i]),]#base::subset(obsData,ID == ids[i])
  }
  
  if(is.null(mov.model)){
    mov.model<-vector('list',length(ids))
    for(i in 1:length(ids)){
      mov.model[[i]] <- ~1
    }
  }
  
  if(is.null(err.model)){
    err.model<-vector('list',length(ids))
  }
  
  if(is.null(activity)){
    activity<-vector('list',length(ids))
  }
  
  if(is.null(drift)){
    drift<-vector('list',length(ids))
    for(i in 1:length(ids)){
      drift[[i]] <- FALSE
    }
  }
  
  if(is.null(constr)){
    constr<-vector('list',length(ids))
    for(i in 1:length(ids)){
      constr[[i]]<-list(lower = -Inf,upper = Inf)
    }
  }
  
  if(is.null(prior)){
    prior<-vector('list',length(ids))
  }
  
  registerDoParallel(cores=ncores) 
  model_fits <- 
    foreach(i = 1:length(ids)) %dopar% {
      
      fit <- crawl::crwMLE(
        mov.model =  mov.model[[i]],
        err.model = err.model[[i]],
        activity = activity[[i]],
        drift = drift[[i]],
        data = ind_data[[i]],
        coord = coord,
        Time.name = Time.name,
        initial.state = initial.state[[i]],
        theta = theta[[i]],
        fixPar = fixPar[[i]],
        method = methodMLE,
        control = control,
        constr = constr[[i]],
        prior = prior[[i]],
        need.hess = need.hess,
        initialSANN = initialSANN,
        attempts = attempts
      )
      
    }
  stopImplicitCluster()
  
  if(is.null(predTime)){
    predTime<-vector('list',length(ids))
    for(i in 1:length(ids)){
      iTime <- range(obsData[which(obsData$ID==ids[i]),][[Time.name]])#range(subset(obsData,ID==i)[[Time.name]])
      predTime[[i]] <- seq(ceiling(iTime[1]),floor(iTime[length(iTime)]),timeStep)
    }
  }
  
  registerDoParallel(cores=ncores)
  predData <- foreach(i = 1:length(ids), .combine = rbind) %dopar% {
    mf<-model_fits[[i]]
    if(!is.null(err.model)) mf$rho[is.na(mf$rho)]<-0
    tmp = crawl::crwPredict(mf, predTime=predTime[[i]])
  }
  stopImplicitCluster()
  
  registerDoParallel(cores=ncores)
  crwSim <- foreach(i = 1:length(model_fits)) %dopar% {
    mf<-model_fits[[i]]
    if(!is.null(err.model)) mf$rho[is.na(mf$rho)]<-0
    tmp = crawl::crwSimulator(mf,predTime=predData[[Time.name]][which(predData$ID==i & predData$locType=="p")], method = methodSim, parIS = parIS,
                              df = df, grid.eps = grid.eps, crit = crit, scale = scale, force.quad)
  }
  stopImplicitCluster()
  
  return(list(model_fits=model_fits,predData=predData,crwSim=crwSim))
}
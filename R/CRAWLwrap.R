#' @export
#' @importFrom sp coordinates
#' @importFrom crawl crwMLE crwPredict
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats formula

CRAWLwrap<-function(obsData, timeStep=1, ncores, 
                    mov.model = NULL, err.model = NULL, activity = NULL, drift = NULL, 
                    coord = c("x", "y"), Time.name = "time", initial.state, theta, fixPar, 
                    method = "L-BFGS-B", control = NULL, constr = NULL, 
                    prior = NULL, need.hess = TRUE, initialSANN = list(maxit = 200), attempts = 1,
                    predTime = NULL)
{
  
  if(!is.character(coord) | length(coord)!=2) stop('coord must be character vector of length 2')
  if(any(!(c("ID",Time.name,coord) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name,coord)[!(c("ID",Time.name,coord) %in% names(obsData))],collapse=","))
  
  toProj = obsData[!is.na(obsData[[coord[1]]]) & !is.na(obsData[[coord[2]]]),c("ID",Time.name,coord)]
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
  
  cat('Fitting',length(ids),'tracks using crawl::crwMLE...')
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
        method = method,
        control = control,
        constr = constr[[i]],
        prior = prior[[i]],
        need.hess = need.hess,
        initialSANN = initialSANN,
        attempts = attempts
      )
      
    }
  stopImplicitCluster()
  cat("DONE\n")
  
  for(i in which(!unlist(lapply(model_fits,function(x) inherits(x,"crwFit"))))){
    warning('crawl::crwMLE for individual ',ids[i],' failed;\n',model_fits[[i]],"   Check crawl::crwMLE arguments and/or consult crawl documentation.")
  }
  
  convFits<-which(unlist(lapply(model_fits,function(x) inherits(x,"crwFit"))))
  model_fits<-model_fits[convFits]
  
  if(is.null(predTime)){
    predTime<-vector('list',length(ids))
    for(i in 1:length(ids)){
      iTime <- range(obsData[which(obsData$ID==ids[i]),][[Time.name]])#range(subset(obsData,ID==i)[[Time.name]])
      predTime[[i]] <- seq(ceiling(iTime[1]),floor(iTime[length(iTime)]),timeStep)
    }
  }
  
  cat('Predicting locations (and uncertainty) for',length(ids),'tracks using crawl::crwPredict...')
  registerDoParallel(cores=ncores)
  predData <- 
    foreach(i = 1:length(convFits), .export="crwPredict", .combine = rbind, .errorhandling="remove") %dopar% {
      crawl::crwPredict(model_fits[[i]], predTime=predTime[[convFits[i]]])
    }
  stopImplicitCluster()
  cat("DONE\n")
  
  return(crwData(list(crwFits=model_fits,crwPredict=predData)))
}
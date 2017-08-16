
#' Fit and predict tracks for using crawl
#'
#' Wrapper function for fitting crawl::crwMLE models and predicting locations with crawl::crwPredict for multiple individuals.
#'
#' @param obsData data.frame object containing fields for animal ID ('ID'), time of observation (identified by \code{Time.name}, must be numeric or POSIXct), 
#' and observed locations (x- and y- coordinates identified by \code{coord}), such as that returned by \code{\link{simData}} when temporally-irregular observed locations or
#' measurement error are included. Alternatively, a 'SpatialPointsDataFrame' object from the package 'sp' will 
#' also be accepted, in which case the \code{coord} values will be taken from the spatial data set and ignored in the arguments.  
#' Note that \code{\link[crawl]{crwMLE}} requires that longitude/latitude coordinates be projected to UTM (i.e., easting/northing). For further details see \code{\link[crawl]{crwMLE}}.
#' @param timeStep Length of the time step at which to predict regular locations from the fitted model. Unless \code{predTime} is specified, the sequence of times
#' is \code{seq(a_i,b_i,timeStep)} where a_i and b_i are the times of the first and last observations for individual i. \code{timeStep} can be numeric (regardless of
#' whether \code{obsData[[Time.name]]} is numeric or POSIXct) or a character string (if \code{obsData[[Time.name]]} is of class POSIXct) containing one of "sec", "min", "hour", "day", "DSTday", "week", "month", "quarter" or "year". 
#' This can optionally be preceded by a positive integer and a space, or followed by "s" (e.g., ``2 hours''; see \code{\link[base]{seq.POSIXt}}). \code{timeStep} is not used for individuals for which \code{predTime} is specified.
#' @param ncores Number of cores to use for parallel processing.
#' @param retryFits Number of times to attempt to achieve convergence and valid (i.e., not NaN) variance estimates after the initial model fit. \code{retryFits} differs
#' from \code{attempts} because \code{retryFits} iteratively uses random perturbations of the current parameter estimates as the initial values for likelihood optimization, while 
#' \code{attempts} uses the same initial values (\code{theta}) for each attempt. 
#' @param mov.model List of mov.model objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one movement model is provided, then the same movement model is used
#' for each individual.
#' @param err.model List of err.model objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one error model is provided, then the same error model is used
#' for each individual (in which case the names of the \code{err.model} components corresponding to easting/longitudinal and northing/latitudinal location error must match \code{coord}).
#' @param activity List of activity objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one activity covariate is provided, then the same activity covariate is used
#' for each individual.
#' @param drift List of drift objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one drift component is provided, then the same drift component is used
#' for each individual.
#' @param coord A 2-vector of character values giving the names of the "x" and
#' "y" coordinates in \code{data}. See \code{\link[crawl]{crwMLE}}.
#' @param Time.name Character indicating name of the location time column.  See \code{\link[crawl]{crwMLE}}.
#' @param initial.state List of initial.state objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one initial state is provided, then the same initial states are used
#' for each individual.
#' @param theta List of theta objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one theta is provided, then the same starting values are used
#' for each individual.
#' @param fixPar List of fixPar objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one fixPar is provided, then the same parameters are held fixed to the given value
#' for each individual.
#' @param method Optimization method that is passed to \code{\link{optim}}.
#' @param control Control list which is passed to \code{\link{optim}}.
#' @param constr List of constr objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one constr is provided, then the same box constraints for the parameters are used
#' for each individual.
#' @param prior List of prior objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one prior is provided, then the same prior is used
#' for each individual.
#' @param need.hess A logical value which decides whether or not to evaluate
#' the Hessian for parameter standard errors
#' @param initialSANN Control list for \code{\link{optim}} when simulated
#' annealing is used for obtaining start values. See details
#' @param attempts The number of times likelihood optimization will be
#' attempted using \code{theta} as the starting values.  Note this is not the same as \code{retryFits}.
#' @param predTime List of predTime objects (see \code{\link[crawl]{crwPredict}}) containing an element for each individual. \code{predTime} can 
#' be specified as an alternative to the automatic sequences generated according to \code{timeStep}.  If only one predTime object is provided, then the same prediction times are used
#' for each individual.
#' @param fillCols Logical indicating whether or not to use the crawl::\code{\link[crawl]{fillCols}} function for filling in missing values in \code{obsData} for which
#' there is a single unique value. Default: FALSE. If the output from \code{crawlWrap} is intended for analyses using \code{\link{fitHMM}} or \code{\link{MIfitHMM}}, 
#' setting \code{fillCols=TRUE} should typically be avoided.
#' 
#' @return A \code{\link{crwData}} object, i.e. a list of:
#' \item{crwFits}{A list of \code{crwFit} objects returned by crawl::crwMLE. See \code{\link[crawl]{crwMLE}}}
#' \item{crwPredict}{A \code{crwPredict} data frame with \code{obsData} merged with the predicted locations. See \code{\link[crawl]{crwPredict}}.}
#' The \code{\link{crwData}} object is used in \code{\link{MIfitHMM}} analyses that account for temporal irregularity or location measurement error.
#' 
#' @details
#' \itemize{
#' \item Consult \code{\link[crawl]{crwMLE}} and \code{\link[crawl]{crwPredict}} for futher details about model fitting and prediction. 
#' 
#' \item Note that the names of the list elements corresponding to each individual in \code{mov.model}, \code{err.model}, \code{activity}, \code{drift}, \code{initial.state},
#' \code{theta}, \code{fixPar}, \code{constr}, \code{prior}, and \code{predTime} must match the individual IDs in \code{obsData}.  If only one element is provided
#' for any of these arguments, then the same element will be applied to all individuals.
#' }
#' 
#' @seealso \code{\link{MIfitHMM}}, \code{\link{simData}}
#' 
#' @examples
#' \dontrun{
#' # extract simulated obsData from example data
#' obsData <- miExample$obsData
#' 
#' # extract crwMLE inputs from example data
#' inits <- miExample$inits # initial state
#' err.model <- miExample$err.model # error ellipse model
#'
#' # Fit crwMLE models to obsData and predict locations 
#' # at default intervals for both individuals
#' crwOut1 <- crawlWrap(obsData=obsData,ncores=1,
#'          theta=c(4,0),fixPar=c(1,1,NA,NA),
#'          initial.state=inits,
#'          err.model=err.model,attempts=100)
#'
#' # Fit the same crwMLE models and predict locations 
#' # at same intervals but specify for each individual using lists
#' crwOut2 <- crawlWrap(obsData=obsData,ncores=1,
#'          theta=list(c(4,0),c(4,0)), fixPar=list(c(1,1,NA,NA),c(1,1,NA,NA)),
#'          initial.state=list(inits,inits),
#'          err.model=list(err.model,err.model),
#'          predTime=list('1'=seq(1,633),'2'=seq(1,686)))
#' }
#' 
#' @export
#' @importFrom crawl crwMLE crwPredict
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats formula
#' @importFrom sp coordinates
#' @importFrom utils capture.output

crawlWrap<-function(obsData, timeStep=1, ncores, retryFits = 0,
                    mov.model = ~1, err.model = NULL, activity = NULL, drift = NULL, 
                    coord = c("x", "y"), Time.name = "time", initial.state, theta, fixPar, 
                    method = "L-BFGS-B", control = NULL, constr = NULL, 
                    prior = NULL, need.hess = TRUE, initialSANN = list(maxit = 200), attempts = 1,
                    predTime = NULL, fillCols = FALSE)
{
  
  if(is.data.frame(obsData)){
    if(!is.character(coord) | length(coord)!=2) stop('coord must be character vector of length 2')
    if(any(!(c("ID",Time.name,coord) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name,coord)[!(c("ID",Time.name,coord) %in% names(obsData))],collapse=","))
  } else if(inherits(obsData,"SpatialPoints")){
    if(any(!(c("ID",Time.name) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name)[!(c("ID",Time.name) %in% names(obsData))],collapse=","))
    coord <- colnames(sp::coordinates(obsData))
  } else stop("obsData must be a data frame or a SpatialPointsDataFrame")
    
  if(retryFits<0) stop("retryFits must be non-negative")
  
  ids = as.character(unique(obsData$ID))
  ind_data<-list()
  
  for(i in ids){
    ind_data[[i]] = obsData[which(obsData$ID==i),]
    if(any(is.na(ind_data[[i]][[Time.name]]))) stop("obsData$",Time.name," cannot contain missing values")
  }
  
  if(is.null(mov.model)){
    mov.model<-list()
    for(i in ids){
      mov.model[[i]] <- ~1
    }
  } else if(is.formula(mov.model)){
    tmpmov.model<-mov.model
    mov.model<-list()
    for(i in ids){
      mov.model[[i]] <- tmpmov.model
    }    
  }
  if(!is.null(names(mov.model))) {
    if(!all(names(mov.model) %in% ids)) stop("mov.model names must match obsData$ID")
    mov.model <- mov.model[ids]
  } else names(mov.model) <- ids
  
  if(is.null(err.model)){
    err.model<-vector('list',length(ids))
    names(err.model)<-ids
  } else if(is.list(err.model)){
    if(all(unlist(lapply(err.model,is.formula)))){
      tmperr.model<-err.model
      err.model<-list()
      for(i in ids){
        err.model[[i]] <- tmperr.model
      } 
    }
  }
  if(!is.null(names(err.model)))  {
    if(!all(names(err.model) %in% ids)) stop("err.model names must match obsData$ID")
    err.model <- err.model[ids]
  } else names(err.model) <- ids
  
  if(is.null(activity)){
    activity<-vector('list',length(ids))
    names(activity)<-ids
  } else if(is.formula(activity)){
    tmpactivity<-activity
    activity<-list()
    for(i in ids){
      activity[[i]] <- tmpactivity
    }    
  }
  if(!is.null(names(activity))) {
    if(!all(names(activity) %in% ids)) stop("activity names must match obsData$ID")
    activity <- activity[ids]
  } else names(activity) <- ids
  
  if(is.null(drift)){
    drift<-list()
    for(i in ids){
      drift[[i]] <- FALSE
    }
  } else if(is.logical(drift)){
    tmpdrift<-drift
    drift<-list()
    for(i in ids){
      drift[[i]] <- tmpdrift
    }    
  }
  if(!is.null(names(drift))) {
    if(!all(names(drift) %in% ids)) stop("drift names must match obsData$ID")
    drift <- drift[ids]
  } else names(drift) <- ids
  
  if(!is.null(names(initial.state))){
    if(all(names(initial.state) %in% c("a","P"))){
      tmpinitial.state<-initial.state
      initial.state<-list()
      for(i in ids){
        initial.state[[i]] <- tmpinitial.state
      } 
    }
  }
  if(!is.null(names(initial.state)))  {
    if(!all(names(initial.state) %in% ids)) stop("initial.state names must match obsData$ID")
    initial.state <- initial.state[ids]
  } else names(initial.state) <- ids
  
  if(!is.list(theta)){
    tmptheta<-theta
    theta<-list()
    for(i in ids){
      theta[[i]] <- tmptheta
    } 
  }
  if(!is.null(names(theta))) {
    if(!all(names(theta) %in% ids)) stop("theta names must match obsData$ID")
    theta <- theta[ids]
  } else names(theta) <- ids
  
  if(!is.list(fixPar)){
    tmpfixPar<-fixPar
    fixPar<-list()
    for(i in ids){
      fixPar[[i]] <- tmpfixPar
    } 
  }
  if(!is.null(names(fixPar))) {
    if(!all(names(fixPar) %in% ids)) stop("fixPar names must match obsData$ID")
    fixPar <- fixPar[ids]
  } else names(fixPar) <- ids
  
  if(is.null(constr)){
    constr<-list()
    for(i in ids){
      constr[[i]]<-list(lower = -Inf,upper = Inf)
    }
  } else if(is.list(constr)){
    if(!is.null(names(constr))){
      if(all(names(constr) %in% c("upper","lower"))){
        tmpconstr<-constr
        constr<-list()
        for(i in ids){
          constr[[i]] <- tmpconstr
        } 
      }
    }
  }
  if(!is.null(names(constr))) {
    if(!all(names(constr) %in% ids)) stop("constr names must match obsData$ID")
    constr <- constr[ids]
  } else names(constr) <- ids
  
  if(is.null(prior)){
    prior<-vector('list',length(ids))
    names(prior)<-ids
  } else if(is.function(prior)){
    tmpprior<-prior
    prior<-list()
    for(i in ids){
      prior[[i]] <- tmpprior
    }    
  }
  if(!is.null(names(prior))) {
    if(!all(names(prior) %in% ids)) stop("prior names must match obsData$ID")
    prior <- prior[ids]
  } else names(prior) <- ids
  
  if(is.null(predTime)){
    predTime<-list()
  }
  if(!is.list(predTime)){
    tmpPredTime<-predTime
    predTime<-list()
    for(i in ids){
      predTime[[i]]<-tmpPredTime
    }
  }
  for(i in ids){
    if(is.null(predTime[[i]])){
      iTime <- range(obsData[which(obsData$ID==i),][[Time.name]])
      predTime[[i]] <- seq(iTime[1],iTime[2],timeStep)
    }
  }
  if(!is.null(names(predTime))) {
    if(!all(names(predTime) %in% ids)) stop("predTime names must be character strings that match elements of obsData$ID. See example in ?crawlWrap.")
    predTime <- predTime[ids]
  } else names(predTime) <- ids
  
  cat('Fitting',length(ids),'track(s) using crawl::crwMLE...')
  registerDoParallel(cores=ncores) 
  model_fits <- 
    foreach(i = ids) %dopar% {
      
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
  
  names(model_fits) <- ids
  
  # Check crwFits and re-try based on retryFits
  for(i in ids){
    if(!inherits(model_fits[[i]],"crwFit"))
      warning('crawl::crwMLE for individual ',i,' failed;\n',model_fits[[i]],"   Check crawl::crwMLE arguments and/or consult crawl documentation.")
    else {
      if((model_fits[[i]]$convergence | any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))]))) | retryFits){
        if(retryFits){
          fitCount<-0
          fit <- model_fits[[i]]
          if(model_fits[[i]]$convergence | any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))]))){
            if(model_fits[[i]]$convergence)
              cat('crawl::crwMLE for individual',i,'has suspect convergence: ',model_fits[[i]]$message,"\n")
            if(any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))])))
              cat('crawl::crwMLE for individual',i,'has NaN variance estimate(s)\n')
            cat('Attempting to achieve convergence and valid variance estimates for individual ',i,". Press 'esc' to force exit from 'crawlWrap'\n",sep="")
          } else {
            cat('Attempting to improve fit for individual ',i,". Press 'esc' to force exit from 'crawlWrap'\n",sep="")
          }
          curFit<-fit
          while(fitCount<retryFits){ # & (fit$convergence | any(is.na(fit$se[which(is.na(fixPar[[i]]))])))){
            cat("\r    Attempt ",fitCount+1," of ",retryFits," -- current log-likelihood value: ",curFit$loglik,"  ...",sep="")
            tmp <- tryCatch(suppressWarnings(suppressMessages(crawl::crwMLE(
              mov.model =  mov.model[[i]],
              err.model = err.model[[i]],
              activity = activity[[i]],
              drift = drift[[i]],
              data = ind_data[[i]],
              coord = coord,
              Time.name = Time.name,
              initial.state = initial.state[[i]],
              theta = fit$estPar + rnorm(length(fit$estPar),0,0.5),
              fixPar = fixPar[[i]],
              method = method,
              control = list(maxit = control$maxit),
              constr = constr[[i]],
              prior = prior[[i]],
              need.hess = need.hess,
              initialSANN = list(maxit = initialSANN$maxit),
              attempts = 1
            ))),error=function(e) e)
            if(inherits(tmp,"crwFit")){
              if(tmp$convergence==0){
                if(tmp$aic < model_fits[[i]]$aic | all(!is.na(tmp$se[which(is.na(fixPar[[i]]))])))
                  fit<-tmp
                if(tmp$aic <= model_fits[[i]]$aic & all(!is.na(tmp$se[which(is.na(fixPar[[i]]))])))
                  curFit<-tmp
              }
            }
            fitCount <- fitCount + 1
          }
          if(curFit$convergence | any(is.na(curFit$se[which(is.na(fixPar[[i]]))]))){
            cat("FAILED\n")
          } else {
            cat("DONE\n")
          }
          model_fits[[i]]<-curFit
        } else {
          if(model_fits[[i]]$convergence)
            warning('crawl::crwMLE for individual ',i,' has suspect convergence: ',model_fits[[i]]$message)
          if(any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))])))
            warning('crawl::crwMLE for individual ',i,' has NaN variance estimate(s)')
        }
      }
    }
  }
  
  convFits <- ids[which(unlist(lapply(model_fits,function(x) inherits(x,"crwFit"))))]
  model_fits <- model_fits[convFits]
  
  if(inherits(obsData[[Time.name]],"POSIXct")){
    td <- utils::capture.output(difftime(predTime[[1]][2],predTime[[1]][1],units="auto"))
    cat('Predicting locations (and uncertainty) at',substr(td,20,nchar(td)),'time steps for',length(convFits),'track(s) using crawl::crwPredict...')
  } else cat('Predicting locations (and uncertainty) for',length(convFits),'track(s) using crawl::crwPredict...')
  
  registerDoParallel(cores=ncores)
  predData <- 
    foreach(i = convFits, .export="crwPredict", .combine = rbind, .errorhandling="remove") %dopar% {
      pD<-crawl::crwPredict(model_fits[[i]], predTime=predTime[[i]])
      pD[[Time.name]][which(pD$locType=="p")]<-predTime[[i]][predTime[[i]]>=min(model_fits[[i]]$data[[Time.name]])]
      if(!fillCols){
        for(j in names(pD)[names(pD) %in% names(ind_data[[i]])]){
          if(!(j %in% c(Time.name,"ID",coord))){
            if(!isTRUE(all.equal(pD[[j]],ind_data[[i]][[j]]))) {
              pD[[j]][pD[[Time.name]] %in% ind_data[[i]][[Time.name]]] <- ind_data[[i]][[j]]
              pD[[j]][!(pD[[Time.name]] %in% ind_data[[i]][[Time.name]])] <- NA
            }
          }
        }
      }
      pD
    }
  stopImplicitCluster()
  cat("DONE\n")
  
  return(crwData(list(crwFits=model_fits,crwPredict=predData)))
}
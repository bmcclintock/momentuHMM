
#' Fit and predict tracks for using crawl
#'
#' Wrapper function for fitting crawl::crwMLE models and predicting locations with crawl::crwPredict for multiple individuals.
#'
#' @param obsData data.frame object containing fields for animal ID ('ID'), time of observation (identified by \code{Time.name}, must be numeric or POSIXct), 
#' and observed locations (x- and y- coordinates identified by \code{coord}), such as that returned by \code{\link{simData}} when temporally-irregular observed locations or
#' measurement error are included. Alternatively, a \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sf]{sf}} object will 
#' also be accepted, in which case the \code{coord} values will be taken from the spatial data set and ignored in the arguments.  
#' Note that \code{\link[crawl]{crwMLE}} requires that longitude/latitude coordinates be projected to UTM (i.e., easting/northing). For further details see \code{\link[crawl]{crwMLE}}.
#' @param timeStep Length of the time step at which to predict regular locations from the fitted model. Unless \code{predTime} is specified, the sequence of times
#' is \code{seq(a_i,b_i,timeStep)} where a_i and b_i are the times of the first and last observations for individual i. \code{timeStep} can be numeric (regardless of
#' whether \code{obsData[[Time.name]]} is numeric or POSIXct) or a character string (if \code{obsData[[Time.name]]} is of class POSIXct) containing one of "sec", "min", "hour", "day", "DSTday", "week", "month", "quarter" or "year". 
#' This can optionally be preceded by a positive integer and a space, or followed by "s" (e.g., ``2 hours''; see \code{\link[base]{seq.POSIXt}}). \code{timeStep} is not used for individuals for which \code{predTime} is specified.
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing).
#' @param retryFits Number of times to attempt to achieve convergence and valid (i.e., not NaN) variance estimates after the initial model fit.
#' @param retrySD An optional list of scalars or vectors for each individual indicating the standard deviation to use for normal perturbations of \code{theta} when \code{retryFits>0} (or \code{attempts>1}). 
#' Instead of a list object, \code{retrySD} can also be a scalar or a vector, in which case the same values are used for each each individual.
#' If a scalar is provided, then the same value is used for each parameter. If a vector is provided, it must be of length \code{length(theta)} for the corresponding individual(s). Default: 1, i.e., a standard deviation of 1 is used
#' for all parameters of all individuals. Ignored unless \code{retryFits>0} (or \code{attempts>1}).
#' @param retryParallel Logical indicating whether or not to perform \code{retryFits} attempts for each individual in parallel. Default: FALSE. Ignored unless \code{retryFits>0} and \code{ncores>1}.
#' Note that when attempts are done in parallel (i.e. \code{retryParallel=TRUE}), the current value for the log-likelihood of each individual and warnings about convergence are not printed to the console.
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
#' @param proj A list of valid epsg integer codes or proj4string for \code{obsData} that does not
#' inherit either 'sf' or 'sp'. A valid 'crs' list is also accepted. Otherwise, ignored. If only one proj is provided, then the same projection is used
#' for each individual.
#' @param Time.name Character indicating name of the location time column.  See \code{\link[crawl]{crwMLE}}.
#' @param time.scale character. Scale for conversion of POSIX time to numeric for modeling. Defaults to "hours".
#' @param theta List of theta objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one theta is provided, then the same starting values are used
#' for each individual. If theta is not specified, then \code{\link[crawl]{crwMLE}} default values are used (i.e. each parameter is started at zero).
#' @param fixPar List of fixPar objects (see \code{\link[crawl]{crwMLE}}) containing an element for each individual. If only one fixPar is provided, then the same parameters are held fixed to the given value
#' for each individual. If fixPar is not specified, then no parameters are fixed.
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
#' attempted in cases where the fit does not converge or is otherwise non-valid. Note this is not the same as \code{retryFits} because \code{attempts} only applies when the current fit clearly does not appear to have converged; \code{retryFits} will proceed with additional model fitting attempts regardless of the model output.
#' @param predTime List of predTime objects (see \code{\link[crawl]{crwPredict}}) containing an element for each individual. \code{predTime} can 
#' be specified as an alternative to the automatic sequences generated according to \code{timeStep}.  If only one predTime object is provided, then the same prediction times are used
#' for each individual.
#' @param fillCols Logical indicating whether or not to use the crawl::\code{\link[crawl]{fillCols}} function for filling in missing values in \code{obsData} for which
#' there is a single unique value. Default: FALSE. If the output from \code{crawlWrap} is intended for analyses using \code{\link{fitHMM}} or \code{\link{MIfitHMM}}, 
#' setting \code{fillCols=TRUE} should typically be avoided.
#' @param coordLevel Character string indicating the level of the hierarchy for the location data. Ignored unless \code{obsData} includes a 'level' field.
#' @param ... Additional arguments that are ignored.
#' 
#' @return A \code{\link{crwData}} or \code{\link{crwHierData}} object, i.e. a list of:
#' \item{crwFits}{A list of \code{crwFit} objects returned by crawl::crwMLE. See \code{\link[crawl]{crwMLE}}}
#' \item{crwPredict}{A \code{crwPredict} data frame with \code{obsData} merged with the predicted locations. See \code{\link[crawl]{crwPredict}}.}
#' The \code{\link{crwData}} object is used in \code{\link{MIfitHMM}} analyses that account for temporal irregularity or location measurement error.
#' 
#' @details
#' \itemize{
#' \item Consult \code{\link[crawl]{crwMLE}} and \code{\link[crawl]{crwPredict}} for futher details about model fitting and prediction. 
#' 
#' \item Note that the names of the list elements corresponding to each individual in \code{mov.model}, \code{err.model}, \code{activity}, \code{drift},
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
#' # error ellipse model
#' err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
#'
#' # Fit crwMLE models to obsData and predict locations 
#' # at default intervals for both individuals
#' crwOut1 <- crawlWrap(obsData=obsData,
#'          theta=c(4,0),fixPar=c(1,1,NA,NA),
#'          err.model=err.model,attempts=100)
#'
#' # Fit the same crwMLE models and predict locations 
#' # at same intervals but specify for each individual using lists
#' crwOut2 <- crawlWrap(obsData=obsData,
#'          theta=list(c(4,0),c(4,0)), fixPar=list(c(1,1,NA,NA),c(1,1,NA,NA)),
#'          err.model=list(err.model,err.model),
#'          predTime=list('1'=seq(1,633),'2'=seq(1,686)))
#' }
#' 
#' @export
#' @importFrom crawl crwMLE crwPredict displayPar
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom lubridate with_tz
#' @importFrom doRNG %dorng%
#' @importFrom stats formula setNames
#' @importFrom sp coordinates
#' @importFrom utils capture.output

crawlWrap<-function(obsData, timeStep=1, ncores = 1, retryFits = 0, retrySD = 1, retryParallel = FALSE,
                    mov.model = ~1, err.model = NULL, activity = NULL, drift = NULL, 
                    coord = c("x", "y"), proj = NULL, Time.name = "time", time.scale = "hours", theta, fixPar, 
                    method = "L-BFGS-B", control = NULL, constr = NULL, 
                    prior = NULL, need.hess = TRUE, initialSANN = list(maxit = 200), attempts = 1,
                    predTime = NULL, fillCols = FALSE, coordLevel = NULL, ...)
{
  
  if(is.data.frame(obsData)){
    if(!inherits(obsData,"sf")) {
      if(!is.character(coord) | length(coord)!=2) stop('coord must be character vector of length 2')
      if(any(!(c("ID",Time.name,coord) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name,coord)[!(c("ID",Time.name,coord) %in% names(obsData))],collapse=","))
      if(any(coord %in% c("mu.x","nu.x","mu.y","nu.y","se.mu.x","se.nu.x","se.mu.y","se.nu.y","speed"))) stop("coordinates cannot include the following names: mu.x, nu.x, mu.y, nu.y, se.mu.x, se.nu.x, se.mu.y, se.nu.y, or speed \n   please choose different coord names")
    } else {
      if(any(!(c("ID",Time.name) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name)[!(c("ID",Time.name) %in% names(obsData))],collapse=","))
    }
  } else if(inherits(obsData,"SpatialPoints")){
    if(any(!(c("ID",Time.name) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name)[!(c("ID",Time.name) %in% names(obsData))],collapse=","))
    coord <- colnames(sp::coordinates(obsData)) # c("x","y") 
  } else stop("obsData must be a data frame or a SpatialPointsDataFrame")
    
  if(retryFits<0) stop("retryFits must be non-negative")
  if(attempts<1) stop("attempts must be >=1")
  
  ids = as.character(unique(obsData$ID))
  ind_data<-list()
  
  hierInd <- FALSE
  for(i in ids){
    ind_data[[i]] = obsData[which(obsData$ID==i),]
    if(any(is.na(ind_data[[i]][[Time.name]]))) stop("obsData$",Time.name," cannot contain missing values")
    #if(any(duplicated(ind_data[[i]][[Time.name]]))){
    #  if(is.null(ind_data[[i]]$level) | is.null(coordLevel)) stop("duplicated times can only be included when coordLevel is specified and obsData includes a 'level' field")
    #  if(!is.factor(ind_data[[i]]$level)) stop("'level' field must be a factor")
    #  if(!(coordLevel %in% levels(ind_data[[i]]$level))) stop("'coordLevel' not found in 'level' field")
    #  ind_data[[i]] <- obsData[which(obsData$ID==i & obsData$level==coordLevel),]
    #  hierInd <- TRUE
    #} else {
      if(!is.null(ind_data[[i]]$level)){
        if(is.factor(ind_data[[i]]$level) & is.null(coordLevel)) stop("'level' field cannot be a factor unless 'coordLevel' is specified")
        if(!is.factor(ind_data[[i]]$level) & !is.null(coordLevel)) stop("'level' field must be a factor when 'coordLevel' is specified")
        if(!is.character(coordLevel) | length(coordLevel)!=1) stop("coordLevel must be a character string")
        if(!(coordLevel %in% levels(ind_data[[i]]$level))) stop("'coordLevel' not found in 'level' field")
        ind_data[[i]] <- obsData[which(obsData$ID==i & obsData$level==coordLevel),]
        hierInd <- TRUE
      } else if(!is.null(coordLevel)) stop("coordLevel can not be specified unless obsData includes a 'level' field")
    #}
  }
  
  nbAnimals <- length(ids)
  
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
    err.model<-vector('list',nbAnimals)
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
    activity<-vector('list',nbAnimals)
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
  
  if(!missing(fixPar)){
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
      for(i in ids){
        if(is.null(fixPar[[i]])){
          fixPar[[i]] <- crawl::displayPar(mov.model =  mov.model[[i]], err.model = err.model[[i]], activity = activity[[i]], drift = drift[[i]], data = ind_data[[i]], Time.name = Time.name)$fixPar
        }
      }
    } else names(fixPar) <- ids
  } else {
    fixPar <- list()
    for(i in ids){
      fixPar[[i]] <- crawl::displayPar(mov.model =  mov.model[[i]], err.model = err.model[[i]], activity = activity[[i]], drift = drift[[i]], data = ind_data[[i]], Time.name = Time.name)$fixPar
    } 
  }
  
  if(!missing(theta)){
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
      for(i in ids){
        if(is.null(theta[[i]])){
          theta[[i]] <- rep(0,sum(is.na(fixPar[[i]])))
        }
      }
    } else names(theta) <- ids
  } else {
    theta <- list()
    for(i in ids){
      theta[[i]] <- rep(0,sum(is.na(fixPar[[i]])))
    }     
  }
  
  if(retryFits>0 | attempts>1){
    if(!is.list(retrySD)){
      tmpretrySD <- retrySD
      retrySD<-list()
      for(i in ids){
        if(length(tmpretrySD)>1){
          if(length(theta[[i]])!=length(tmpretrySD)) stop("retrySD is not the correct length for individual ",i)
          retrySD[[i]] <- tmpretrySD
        } else {
          retrySD[[i]] <- rep(tmpretrySD,length(theta[[i]]))
        }
      }
    } else {
      if(!is.null(names(retrySD))) {
        if(!all(names(retrySD) %in% ids)) stop("retrySD names must match obsData$ID")
        retrySD <- retrySD[ids]
      } else {
        if(length(retrySD)!=nbAnimals) stop('when no list object names are provided, retrySD must be a list of length ',nbAnimals)
        names(retrySD) <- ids
      }
      for(i in ids){
        if(is.null(retrySD[[i]])){
          retrySD[[i]] <- rep(1,length(theta[[i]]))
        } else {
          tmpretrySD <- retrySD[[i]]
          if(length(tmpretrySD)>1){
            if(length(theta[[i]])!=length(tmpretrySD)) stop("retrySD is not the correct length for individual ",i)
          } else {
            retrySD[[i]] <- rep(tmpretrySD,length(theta[[i]]))
          }
        }
      }
    }
    retrySD <- retrySD[ids]
  } else {
    retrySD <- vector('list',nbAnimals)
    names(retrySD) <- ids
    for(i in ids){
      retrySD[[i]] <- rep(1,length(theta[[i]]))
    }
  }
  
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
    prior<-vector('list',nbAnimals)
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
  
  if(is.null(proj)){
    proj<-vector('list',nbAnimals)
    names(proj)<-ids
  } else if(!is.list(proj)){
    tmpproj<-proj
    proj<-list()
    for(i in ids){
      proj[[i]] <- tmpproj
    }    
  }
  if(!is.null(names(proj))) {
    if(!all(names(proj) %in% ids)) stop("proj names must match obsData$ID")
    proj <- proj[ids]
  } else names(proj) <- ids
  
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
      iTime <- range(ind_data[[i]][[Time.name]])
      if(inherits(obsData[[Time.name]],"POSIXct")){
        tzone <- attributes(obsData[[Time.name]])$tzone
        predTime[[i]] <- as.POSIXct(seq(iTime[1],iTime[2],timeStep),tz=tzone)
      } else {
        tzone <- NULL
        predTime[[i]] <- seq(iTime[1],iTime[2],timeStep)
      }
      attributes(predTime[[i]])$tzone <- tzone
    }
  }
  if(!is.null(names(predTime))) {
    if(!all(names(predTime) %in% ids)) stop("predTime names must be character strings that match elements of obsData$ID. See example in ?crawlWrap.")
    predTime <- predTime[ids]
  } else names(predTime) <- ids
  
  cat('Fitting',nbAnimals,'track(s) using crawl::crwMLE...',ifelse(nbAnimals>1 & ncores>1,"","\n"))
  registerDoParallel(cores=ncores)
  withCallingHandlers(model_fits <- 
    foreach(i = ids, .export="crwMLE", .errorhandling="pass", .final = function(x) stats::setNames(x, ids)) %dorng% {
      cat("Individual ",i,"...\n",sep="")
      fit <- crawl::crwMLE(
        data = ind_data[[i]],
        mov.model =  mov.model[[i]],
        err.model = err.model[[i]],
        activity = activity[[i]],
        drift = drift[[i]],
        coord = coord,
        proj = proj[[i]],
        Time.name = Time.name,
        time.scale = time.scale,
        theta = theta[[i]],
        fixPar = fixPar[[i]],
        method = method,
        control = control,
        constr = constr[[i]],
        prior = prior[[i]],
        need.hess = need.hess,
        initialSANN = initialSANN,
        attempts = attempts, 
        retrySD = retrySD[[i]],
        ... = ...
      )
    }
  ,warning=muffleRNGwarning)
  stopImplicitCluster()
  
  convFits <- ids[which(unlist(lapply(model_fits,function(x) inherits(x,"crwFit"))))]
  if(!length(convFits)) stop("crawl::crwMLE failed for all individuals.  Check crawl::crwMLE arguments and/or consult crawl documentation.\n")
    
  cat("DONE\n")
  if(retryFits>0) cat("\n\n")
  
  if(retryParallel & ncores>1){
    for(i in ids){
      if(!inherits(model_fits[[i]],"crwFit"))
        warning('crawl::crwMLE for individual ',i,' failed;\n',model_fits[[i]],"   Check crawl::crwMLE arguments and/or consult crawl documentation.")
    }
  }
  
  # Check crwFits and re-try based on retryFits
  tmpcores <- ncores
  if(!retryParallel) tmpcores <- 1
  else if(ncores>1 & nbAnimals>1) cat("Attempting to achieve convergence and valid variance estimates for each individual in parallel.\n    Press 'esc' to force exit from 'crawlWrap'... ",sep="")
  registerDoParallel(cores=tmpcores)
  withCallingHandlers(model_fits <- foreach(i = ids, .export=c("quietCrawl","crwMLE"), .errorhandling="pass", .final = function(x) stats::setNames(x, ids)) %dorng% {
    if(inherits(model_fits[[i]],"crwFit")){
      if((model_fits[[i]]$convergence | any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))]))) | retryFits){
        if(retryFits){
          fitCount<-0
          fit <- model_fits[[i]]
          if(model_fits[[i]]$convergence | any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))]))){
            if(model_fits[[i]]$convergence)
              cat('\ncrawl::crwMLE for individual',i,'has suspect convergence: ',model_fits[[i]]$message,"\n")
            if(any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))])))
              cat('\ncrawl::crwMLE for individual',i,'has NaN variance estimate(s)\n')
            cat('Attempting to achieve convergence and valid variance estimates for individual ',i,". Press 'esc' to force exit from 'crawlWrap'\n",sep="")
          } else {
            cat('Attempting to improve fit for individual ',i,". Press 'esc' to force exit from 'crawlWrap'\n",sep="")
          }
          curFit<-fit
          while(fitCount<retryFits){ # & (fit$convergence | any(is.na(fit$se[which(is.na(fixPar[[i]]))])))){
            cat("\r    Attempt ",fitCount+1," of ",retryFits," -- current log-likelihood value: ",curFit$loglik,"  ...",sep="")
            tmpFun <- function(){
              tryCatch(suppressWarnings(suppressMessages(
                crawl::crwMLE(data = ind_data[[i]],
                              mov.model =  mov.model[[i]],
                              err.model = err.model[[i]],
                              activity = activity[[i]],
                              drift = drift[[i]],
                              coord = coord,
                              proj = proj[[i]],
                              Time.name = Time.name,
                              time.scale = time.scale,
                              theta = fit$estPar + rnorm(length(fit$estPar),0,retrySD[[i]]),
                              fixPar = fixPar[[i]],
                              method = method,
                              control = control,
                              constr = constr[[i]],
                              prior = prior[[i]],
                              need.hess = need.hess,
                              initialSANN = initialSANN,
                              attempts = 1, 
                              retrySD = retrySD[[i]],
                              ... = ...))),error=function(e){e})}
            tmp <- NULL
            if(retryParallel){
              tmp <- tmpFun()
            } else {
              # hack to suppress crawl cat output
              tmp <- quietCrawl(tmpFun())
            }
            if(inherits(tmp,"crwFit")){
              if(tmp$convergence==0){
                if(tmp$loglik > curFit$loglik | all(!is.na(tmp$se[which(is.na(fixPar[[i]]))])))
                  fit<-tmp
                if(tmp$loglik >= curFit$loglik & all(!is.na(tmp$se[which(is.na(fixPar[[i]]))])))
                  curFit<-tmp
              }
            }
            fitCount <- fitCount + 1
          }
          if(curFit$convergence | any(is.na(curFit$se[which(is.na(fixPar[[i]]))]))){
            message("FAILED\n")
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
    } else {
      warning('\ncrawl::crwMLE for individual ',i,' failed;\n',model_fits[[i]],"   Check crawl::crwMLE arguments and/or consult crawl documentation.\n")
    }
    model_fits[[i]]
  },warning=muffleRNGwarning)
  stopImplicitCluster()
  if(retryParallel & ncores>1 & nbAnimals>1) cat("DONE\n")

  convFits <- ids[which(unlist(lapply(model_fits,function(x) inherits(x,"crwFit"))))]
  if(!length(convFits)) stop("crawl::crwMLE failed for all individuals.  Check crawl::crwMLE arguments and/or consult crawl documentation.\n")
  model_fits <- model_fits[convFits]
  
  txt <- NULL
  if(inherits(obsData[[Time.name]],"POSIXct")){
    td <- list()
    for(i in convFits){
      td[[i]] <- predTime[[i]]
      if(length(predTime[[i]])>1){
        td[[i]] <- utils::capture.output(difftime(predTime[[i]][2],predTime[[i]][1],units="auto"))
        td[[i]] <- substr(td[[i]],20,nchar(td[[i]]))
      }
    }
    if(length(unique(td))==1) txt <- paste('at',td[[1]],'time steps')
  }
  if(retryFits>0) cat("\n")
  cat('\nPredicting locations (and uncertainty)',txt,'for',length(convFits),'track(s) using crawl::crwPredict... ')
  
  if (time.scale %in% c("hours", "hour")) {
    ts <- 60 * 60
  } else if (time.scale %in% c("days", "day")) {
    ts <- 60 * 60 * 24
  } else if (time.scale %in% c("sec", "secs", "second","seconds")) {
    ts <- 1
  } else if (time.scale %in% c("min", "mins", "minute", "minutes")) {
    ts <- 60
  }
  
  registerDoParallel(cores=ncores)
  withCallingHandlers(predData <- 
    foreach(i = convFits, .export="crwPredict", .combine = rbind, .errorhandling="remove") %dorng% {
      pD<-crawl::crwPredict(model_fits[[i]], predTime=predTime[[i]],return.type = "flat")
      if(inherits(ind_data[[i]][[Time.name]],"POSIXct") && attributes(pD[[Time.name]])$tzone != attributes(ind_data[[i]][[Time.name]])$tzone){
        pD[[Time.name]] <- lubridate::with_tz(pD[[Time.name]],tz=attributes(ind_data[[i]][[Time.name]])$tzone)
      }
      if(length(predTime[[i]])>1){
        pD[[Time.name]][which(pD$locType=="p")]<-predTime[[i]][predTime[[i]]>=min(model_fits[[i]]$data[[Time.name]])]
      } else if(inherits(obsData[[Time.name]],"POSIXct")){
        pD[[Time.name]] <- as.POSIXct(pD$TimeNum*ts,origin="1970-01-01 00:00:00",tz=attributes(pD[[Time.name]])$tzone)
      }
      if(!fillCols){
        # remove duplicated observation times (because this what crwPredict does)
        dups <- duplicated(ind_data[[i]][[Time.name]])
        tmpind_data <- as.data.frame(ind_data[[i]][!dups,,drop=FALSE])
        for(j in names(pD)[names(pD) %in% names(ind_data[[i]])]){
          if(!(j %in% c(Time.name,"ID",coord))){
            if(!isTRUE(all.equal(pD[[j]],tmpind_data[[j]]))) {
              pD[[j]][pD[[Time.name]] %in% tmpind_data[[Time.name]]] <- tmpind_data[[j]]
              pD[[j]][!(pD[[Time.name]] %in% tmpind_data[[Time.name]])] <- NA
            }
          }
        }
      }
      if(!is.null(coordLevel)) pD$level <- coordLevel
      pD
    }
  ,warning=muffleRNGwarning)
  stopImplicitCluster()
  
  if(hierInd){
    
    pData <- predData
    predData <- data.frame()
    
    for(i in convFits){
      
      ipData <- pData[which(pData$ID==i),]
      ipData$level <- factor(ipData$level,levels=levels(obsData$level))
      tmpData <- obsData[which(obsData$ID==i),]
      #tmpData[[Time.name]] <- as.POSIXct(as.numeric(tmpData[[Time.name]]),origin="1970-01-01 00:00:00",tz=attributes(pData[[Time.name]])$tzone)
      tmpData <- merge(tmpData,ipData[,c(Time.name,"level"),drop=FALSE],all=TRUE,by=c(Time.name,"level"))
      tmpData$level[is.na(tmpData$level)] <- coordLevel
      tmpData$ID[is.na(tmpData$ID)] <- i
      
      for(jj in names(tmpData)[!(names(tmpData) %in% names(ipData))]){
        ipData[[jj]] <- rep(NA,nrow(ipData))
      }
      for(jj in names(ipData)[!(names(ipData) %in% names(tmpData))]){
        tmpData[[jj]] <- rep(NA,nrow(tmpData))
      }
      
      ipData <- ipData[,names(tmpData)]
      
      #for(jj in 1:nrow(ipData)){
      #    tmpInd <- which(tmpData$time==ipData$time[jj] & tmpData$level==coordLevel)
      #    tmpData[tmpInd,is.na(tmpData[tmpInd,])] <- ipData[jj,is.na(tmpData[tmpInd,])]
      #}
      tmpInd <- which(tmpData$time %in% ipData$time & tmpData$level==coordLevel)
      tmpData[tmpInd,] <- Map(function(x,y) {x[is.na(x)] <- y[is.na(x)]; x}, tmpData[tmpInd,], ipData[ipData$time %in% tmpData$time,])
      
      predData<-rbind(predData,tmpData)
    }
    predData <- predData[,names(pData)]
    attrNames <- names(attributes(pData))[!(names(attributes(pData)) %in% names(attributes(predData)))]
    attributes(predData)[attrNames] <- attributes(pData)[attrNames]
    attr(predData,"coordLevel") <- coordLevel
    class(predData) <- append(c("crwPredict","hierarchical"),class(predData))
    cat("DONE\n")
    return(crwHierData(list(crwFits=model_fits,crwPredict=predData)))    
  } else {
    cat("DONE\n")
    return(crwData(list(crwFits=model_fits,crwPredict=predData)))    
  }
}

quietCrawl <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x))
} 
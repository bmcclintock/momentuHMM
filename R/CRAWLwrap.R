#' @export
#' @importFrom crawl crwMLE crwPredict
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats formula

CRAWLwrap<-function(obsData, timeStep=1, ncores, retryFits = 0,
                    mov.model = ~1, err.model = NULL, activity = NULL, drift = NULL, 
                    coord = c("x", "y"), Time.name = "time", initial.state, theta, fixPar, 
                    method = "L-BFGS-B", control = NULL, constr = NULL, 
                    prior = NULL, need.hess = TRUE, initialSANN = list(maxit = 200), attempts = 1,
                    predTime = NULL)
{
  
  if(!is.character(coord) | length(coord)!=2) stop('coord must be character vector of length 2')
  if(any(!(c("ID",Time.name,coord) %in% names(obsData)))) stop('obsData is missing ',paste(c("ID",Time.name,coord)[!(c("ID",Time.name,coord) %in% names(obsData))],collapse=","))
  if(retryFits<0) stop("retryFits must be non-negative")
  
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
  } else if(is.formula(mov.model)){
    tmpmov.model<-mov.model
    mov.model<-vector('list',length(ids))
    for(i in 1:length(ids)){
      mov.model[[i]] <- tmpmov.model
    }    
  }
  if(!is.null(names(mov.model))) mov.model <- mov.model[ids]
  
  if(is.null(err.model)){
    err.model<-vector('list',length(ids))
  } else if(is.list(err.model)){
    if(!is.null(names(err.model))){
      if(all(names(err.model) %in% c(coord,"rho"))){
        tmperr.model<-err.model
        err.model<-vector('list',length(ids))
        for(i in 1:length(ids)){
          err.model[[i]] <- tmperr.model
        } 
      }
    }
  }
  if(!is.null(names(err.model)))  err.model <- err.model[ids]
  
  if(is.null(activity)){
    activity<-vector('list',length(ids))
  } else if(is.formula(activity)){
    tmpactivity<-activity
    activity<-vector('list',length(ids))
    for(i in 1:length(ids)){
      activity[[i]] <- tmpactivity
    }    
  }
  if(!is.null(names(activity))) activity <- activity[ids]
  
  if(is.null(drift)){
    drift<-vector('list',length(ids))
    for(i in 1:length(ids)){
      drift[[i]] <- FALSE
    }
  } else if(is.logical(drift)){
    tmpdrift<-drift
    drift<-vector('list',length(ids))
    for(i in 1:length(ids)){
      drift[[i]] <- tmpdrift
    }    
  }
  if(!is.null(names(drift))) drift <- drift[ids]
  
  if(!is.list(theta)){
    tmptheta<-theta
    theta<-vector('list',length(ids))
    for(i in 1:length(ids)){
      theta[[i]] <- tmptheta
    } 
  }
  if(!is.null(names(theta))) theta <- theta[ids]
  
  if(!is.list(fixPar)){
    tmpfixPar<-fixPar
    fixPar<-vector('list',length(ids))
    for(i in 1:length(ids)){
      fixPar[[i]] <- tmpfixPar
    } 
  }
  if(!is.null(names(fixPar))) fixPar <- fixPar[ids]
  
  if(is.null(constr)){
    constr<-vector('list',length(ids))
    for(i in 1:length(ids)){
      constr[[i]]<-list(lower = -Inf,upper = Inf)
    }
  } else if(is.list(constr)){
    if(!is.null(names(constr))){
      if(all(names(constr) %in% c("upper","lower"))){
        tmpconstr<-constr
        constr<-vector('list',length(ids))
        for(i in 1:length(ids)){
          constr[[i]] <- tmpconstr
        } 
      }
    }
  }
  if(!is.null(names(constr))) constr <- constr[ids]
  
  if(is.null(prior)){
    prior<-vector('list',length(ids))
  } else if(is.function(prior)){
    tmpprior<-prior
    prior<-vector('list',length(ids))
    for(i in 1:length(ids)){
      prior[[i]] <- tmpprior
    }    
  }
  if(!is.null(names(prior))) prior <- prior[ids]
  
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
  
  # Check crwFits and re-try based on retryFits
  for(i in 1:length(ids)){
    if(!inherits(model_fits[[i]],"crwFit"))
      warning('crawl::crwMLE for individual ',ids[i],' failed;\n',model_fits[[i]],"   Check crawl::crwMLE arguments and/or consult crawl documentation.")
    else {
      if(model_fits[[i]]$convergence | any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))]))){
        if(retryFits){
          fitCount<-0
          fit <- model_fits[[i]]
          if(model_fits[[i]]$convergence)
            cat('crawl::crwMLE for individual',ids[i],'has suspect convergence: ',model_fits[[i]]$message,"\n")
          if(any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))])))
            cat('crawl::crwMLE for individual',ids[i],'has NaN variance estimate(s)\n')
          cat('Attempting to achieve convergence and valid variance estimates for individual ',ids[i],". Press 'esc' to force exit from 'CRAWLwrap'\n",sep="")
          while(fitCount<retryFits & (fit$convergence | any(is.na(fit$se[which(is.na(fixPar[[i]]))])))){
            cat("\r    Attempt ",fitCount+1," of ",retryFits,"...",sep="")
            tmp <- suppressWarnings(suppressMessages(crawl::crwMLE(
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
            )))
            if(inherits(tmp,"crwFit"))
              if(tmp$convergence==0)
                if(tmp$aic < model_fits[[i]]$aic | all(!is.na(tmp$se[which(is.na(fixPar[[i]]))])))
                  fit<-tmp
            fitCount <- fitCount + 1
          }
          cat("DONE\n")
          model_fits[[i]]<-fit
        } else {
          if(model_fits[[i]]$convergence)
            warning('crawl::crwMLE for individual ',ids[i],' has suspect convergence: ',model_fits[[i]]$message)
          if(any(is.na(model_fits[[i]]$se[which(is.na(fixPar[[i]]))])))
            warning('crawl::crwMLE for individual ',ids[i],' has NaN variance estimate(s)')
        }
      }
    }
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
  
  cat('Predicting locations (and uncertainty) for',length(convFits),'track(s) using crawl::crwPredict...')
  registerDoParallel(cores=ncores)
  predData <- 
    foreach(i = 1:length(convFits), .export="crwPredict", .combine = rbind, .errorhandling="remove") %dopar% {
      crawl::crwPredict(model_fits[[i]], predTime=predTime[[convFits[i]]])
    }
  stopImplicitCluster()
  cat("DONE\n")
  
  return(crwData(list(crwFits=model_fits,crwPredict=predData)))
}
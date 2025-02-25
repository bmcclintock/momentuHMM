## code modified from the hmmTMB package: Michelot T, Glennie R (2023). hmmTMB: Fit Hidden Markov Models using Template Model Builder. R package version 1.0.2, <https://CRAN.R-project.org/package=hmmTMB>.

#' @rawNamespace useDynLib(momentuHMM, .registration=TRUE); useDynLib(momentuHMM_TMBExports)
# #' @importFrom methods as
#' @importFrom stats update
# #' @importFrom mgcv gam
# #' @importFrom optimx optimx
fitTMB <- function(data,dist,nbStates,p,estAngleMean,oneInflation,zeroInflation,DM,DMinputs,formula,fixParIndex,stationary,prior,knownStates,betaCons,control,formulaDelta,covsDelta,workBounds,fit,CT,dtIndex,kappa,hessian,crwST=FALSE) {
  
  for(pkg in c("Matrix","optimx","mgcv")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop('Package \"',pkg,'\" needed for optMethod=\"TMB\". Please install it.',
           call. = FALSE)
    }
  }
  
  if(nbStates<2){
    #stop("nbStates must be >1 when optMethod='TMB'")
    fixParIndex$beta0$beta <- matrix(0,1,1)
    fixParIndex$delta0 <- 0
    fixParIndex$fixPar$beta <- NA
    fixParIndex$fixPar$delta <- NA
    workBounds$beta <- matrix(c(-Inf,Inf),1,2)
    workBounds$delta <- matrix(c(-Inf,Inf),1,2)
  }
  nbAnimals <- length(unique(data$ID))
  distnames <- names(dist)
  distcode <- as.vector(sapply(distnames,function(x) dist_code(dist[[x]],zeroInflation[[x]],oneInflation[[x]])))#as.vector(sapply(self$obs()$dists(), function(d) d$code()))
  names(distcode) <- distnames
  
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  
  if(CT){
    dt <- data$dt
  } else dt <- rep(1,nrow(data))
  
  oldpriors <- list(coeff_fe_obs=matrix(NA,length(unlist(fixParIndex$Par0[distnames])),2,dimnames = list(NULL,c("mean","sd"))),
                 coeff_fe_hid=matrix(NA,length(fixParIndex$beta0$beta),2,dimnames = list(NULL,c("mean","sd"))),
                 log_delta0=matrix(NA,length(fixParIndex$delta0),2,dimnames=list(NULL,c("mean","sd"))),
                 log_lambda_obs=matrix(NA,0,2,dimnames=list(NULL,c("mean","sd"))),
                 log_lambda_hid=matrix(NA,0,2,dimnames=list(NULL,c("mean","sd"))))#self$oldpriors()
  
  if(!is.null(prior)){
    if(!is.list(prior) || !length(names(prior))) stop("prior must be a named list")
    parCount <- 0
    for(i in distnames){
      if(!is.null(prior[[i]])){
        if(!is.matrix(prior[[i]]) || !all(dim(prior[[i]])==c(length(fixParIndex$Par0[[i]]),2))) stop("prior$",i," must be a ",length(fixParIndex$Par0[[i]])," x 2 matrix")
        oldpriors$coeff_fe_obs[parCount+1:length(fixParIndex$Par0[[i]]),] <- prior[[i]]
      }
      parCount <- parCount + length(fixParIndex$Par0[[i]])
    }
    if(!is.null(prior$beta)){
      if(!is.matrix(prior$beta) || !all(dim(prior$beta)==c(length(fixParIndex$beta0$beta),2))) stop("prior$beta must be a ",length(fixParIndex$beta0$beta)," x 2 matrix")
      oldpriors$coeff_fe_hid <- prior$beta
    }
  }
  
  Par0 <- fixParIndex$Par0
  DMnames <- list()
  for(i in distnames){
    if(is.null(DM[[i]])) Par0[[i]] <- n2w(fixParIndex$Par0[i],p$bounds,NULL,NULL,nbStates,estAngleMean,DM,p$Bndind,dist,TMB=TRUE)
    DMnames[[i]] <- colnames(DMinputs$fullDM[[i]])
    if(dist[[i]] %in% angledists){
      if(!estAngleMean[[i]]) {
        Par0[[i]] <- c(rep(0,nbStates),Par0[[i]])
        p$parSize[[i]] <- p$parSize[[i]] + 1
        p$parNames[[i]] <- c("mean",p$parNames[[i]])
        fixParIndex$fixPar[[i]] <- c(rep(NA,nbStates),fixParIndex$fixPar[[i]])
        workBounds[[i]] <- rbind(matrix(c(-Inf,Inf),nbStates,2,byrow=TRUE),workBounds[[i]])
        if(is.formula(DM[[i]][[p$parNames[[i]]]])) {
          DM[[i]][["mean"]] <- ~1
          DM[[i]] <- DM[[i]][p$parNames[[i]]]
        }
        DMnames[[i]] <- c(paste0("mean_",1:nbStates,":(Intercept)"),colnames(DMinputs$fullDM[[i]]))
      }
    }
    #if(is.null(DM[[i]]) || is.null(names(Par0[[i]]))) names(Par0[[i]]) <- colnames(DMinputs$fullDM[[i]])
    names(Par0[[i]]) <- DMnames[[i]]
    parNames <- names(DM[[i]])
    if((dist[[i]] %in% rwdists) && dist[[i]]!="ctcrw"){
      dimDM <- length(DM[[i]])
      labs <- lapply(DM[[i]][which(grepl("mean",parNames))],function(x) attr(stats::terms(x,keep.order=TRUE),"term.labels"))
      if(any(grepl("langevin",unlist(DM[[i]][which(grepl("mean",parNames))])))){
        distcode[i] <- 26
      } else {
        if(length(unlist(labs))) distcode[i] <- 27
        else distcode[i] <- 28
      }
      tm1Ind <- tm1AllInd <- NULL
      for(j in which(grepl("mean",parNames))){
        m_tm1 <- paste0(i,gsub("mean","",parNames[j]),"_tm1")
        tm1AllNames <- labs[[j]][which(grepl(m_tm1,labs[[j]]))]
        if(!m_tm1 %in% labs[[j]]) tm1AllNames <- c(m_tm1,tm1AllNames)
        labsj <- labs[[j]][which(!grepl(m_tm1,labs[[j]]))]
        if(any(grepl("langevin",unlist(DM[[i]])))){
          if(!all(grepl("langevin",labsj))) stop("non-Langevin terms included in DM$",i,"$",parNames[j])
        }
        if(distcode[i] %in% c(26,27)){
          tm1Names <- paste0(parNames[j],"_",1:nbStates,":",i,gsub("mean","",parNames[j]),"_tm1")
          tm1Ind <- c(tm1Ind,unlist(lapply(tm1Names,function(x) which(grepl(x,DMnames[[i]])))))
          tm1AllInd <- c(tm1AllInd,unlist(lapply(tm1AllNames,function(x) which(grepl(x,DMnames[[i]])))))
          DM[[i]][[paste0(parNames[j],"_tm1")]] <- stats::as.formula(paste0("~0+",m_tm1,paste0(ifelse(tm1AllNames!=m_tm1,paste0("+",tm1AllNames),""),collapse="")))
          DM[[i]][[parNames[j]]] <- stats::as.formula(paste0("~0",ifelse(length(labsj),paste0("+",paste0(labsj,collapse="+")),"")))
          # reorder Par0
          if(fit && (!all(is.na(fixParIndex$fixPar[[i]][tm1Ind])) | !all(Par0[[i]][c(tm1Ind)]==1))) stop("Parameters corresponding to ",i,".x_tm1 and ",i,".y_tm1 must be fixed to 1")
          Par0[[i]] <- Par0[[i]][c(tm1AllInd,(1:length(DMnames[[i]]))[-tm1AllInd])]
          fixParIndex$fixPar[[i]] <- fixParIndex$fixPar[[i]][c(tm1AllInd,(1:length(DMnames[[i]]))[-tm1AllInd])]
          workBounds[[i]] <- workBounds[[i]][c(tm1AllInd,(1:length(DMnames[[i]]))[-tm1AllInd]),]
          DM[[i]] <- DM[[i]][c(which(grepl("_tm1$",names(DM[[i]]))),which(!grepl("_tm1$",names(DM[[i]]))))]
          p$parNames[[i]] <- names(DM[[i]])
          p$parSize[[i]] <- length(names(DM[[i]]))
          tm1Ind <- match(tm1Ind,tm1AllInd)
          tm1AllInd <- 1:length(tm1AllInd)
        } else {
          DM[[i]][[parNames[j]]] <- stats::as.formula(paste0("~0+",paste0(i,gsub("mean","",parNames[j]),"_tm1")))
          if(fit && (!all(is.na(fixParIndex$fixPar[[i]][tm1Ind])) | !all(Par0[[i]][c(tm1Ind)]==1))) stop("Parameters corresponding to ",i,".x_tm1 and ",i,".y_tm1 (and ",i,".z_tm1 if trivariate normal) must be fixed to 1")
        }
      }
      fixParIndex$fixPar[[i]] <- reorderFix(fixParIndex$fixPar[[i]])
    }
    lag <- 0
    if(is.list(DM[[i]]) | is.null(DM[[i]])) {
      DM[i] <- make_formulas(DM[i],i,p$parNames[i],nbStates)
      if((dist[[i]] %in% rwdists) && (dist[[i]]!="ctcrw")){
        missingPar <- NULL
        for(j in names(DM[[i]])[which(grepl("mean",names(DM[[i]])))]){
          for(s in names(DM[[i]][[j]])){
            formTerms <- stats::terms(DM[[i]][[j]][[s]],keep.order=TRUE)
            lag <- get_crwlag(DM[[i]][[j]][[s]],lag)
            if(lag){
              for(jj in all.vars(DM[[i]][[j]][[s]],data)){
                attr(data[[jj]],"aInd") <- aInd
              }
            }
            if(!length(attr(formTerms,"term.labels")) && !attr(formTerms,"intercept")){
              DM[[i]][[j]][[s]] <- ~1
              #DM[[i]][[paste0(j,"_tm1")]][[s]] <- stats::as.formula(paste0("~0+",i,gsub("mean","",j),"_tm1"))
              missingPar <- c(missingPar,paste0(c(i,j,s),collapse="."))
            } else if(grepl("_tm1",j)){
              tInd <- which(!attr(formTerms,"term.labels") %in% paste0(i,gsub("mean","",j)))
              if(length(tInd)){
                tmpTerms <- gsub(":","*",attr(formTerms,"term.labels")[tInd])
                DM[[i]][[j]][[s]] <- stats::as.formula(paste0("~0+",paste0(i,gsub("mean","",j)),paste0("+I(",tmpTerms," * dt)",collapse="")))
              }
            }
          }
        }
        tmpX_fe <- colnames(make_matrices(DM[i],data)$X_fe)#self$hid()$terms()
        if(length(tmpX_fe)!=length(Par0[[i]])){
          mPar <- mapply(function(x) which(grepl(x,tmpX_fe)),missingPar)
          tmpX_fe[-mPar] <- names(Par0[[i]])
          tmpPar <- tmpFix <- numeric(length(tmpX_fe))
          tmpWB <- matrix(nrow=length(tmpX_fe),ncol=2)
          if(!is.null(prior[[i]])) {
            tmpPrior <- matrix(NA,nrow=length(tmpX_fe),ncol=2)
            rownames(tmpPrior) <- tmpX_fe
          }
          names(tmpPar) <- names(tmpFix) <- tmpX_fe
          rownames(tmpWB) <- tmpX_fe
          tmpPar[match(names(Par0[[i]]),tmpX_fe)] <- Par0[[i]]
          tmpWB[match(names(Par0[[i]]),tmpX_fe),] <- workBounds[[i]] 
          tmpWB[mPar,] <- matrix(c(-Inf,Inf),nrow=length(mPar),ncol=2,byrow=TRUE)
          if(!is.null(prior[[i]])) {
            tmpPrior[match(names(Par0[[i]]),tmpX_fe),] <- prior[[i]]
            prior[[i]] <- tmpPrior
          }
          tmpFix[match(names(Par0[[i]]),tmpX_fe)] <- fixParIndex$fixPar[[i]]
          tmpFix[mPar] <- NA
          Par0[[i]] <- tmpPar
          workBounds[[i]] <- tmpWB
          fixParIndex$fixPar[[i]] <- tmpFix
        }
      }
    } else stop("Pseudo-design matrices are not permitted when optMethod='TMB'. DM$",i," must be specified using formulas")
  }
  distpar <- unlist(p$parSize)#as.vector(sapply(self$obs()$dists(), function(d) d$npar()))
  
  mod_mat_obs <- make_matrices(DM,data)#self$hid()$terms()
  
  X_fe_obs <- mod_mat_obs$X_fe
  X_re_obs <- mod_mat_obs$X_re
  S_obs <- mod_mat_obs$S
  ncol_re_obs <- mod_mat_obs$ncol_re
  
  mod_mat_hid <- getTPMmat(formula,nbStates,fixParIndex$betaRef,data)#self$hid()$terms()
  if(nbStates>=2) {
    X_fe_hid <- mod_mat_hid$X_fe
  } else {
    X_fe_hid <- suppressMessages(methods::as(matrix(1,nrow(data),1), "dgCMatrix"))
    colnames(X_fe_hid) <- "S1>S2.(Intercept)"
  }
  X_re_hid <- mod_mat_hid$X_re
  S_hid <- mod_mat_hid$S
  ncol_re_hid <- mod_mat_hid$ncol_re
  
  ldelta <- get_tmbDelta(fixParIndex$delta0,ifelse(nbStates>1,nbStates,2),nbAnimals,stationary,formulaDelta,covsDelta,fixParIndex$fixPar$delta,workBounds$delta)
  ldelta0 <- ldelta$ldelta0
  workBounds$delta <- ldelta$workBounds
  fixParIndex$delta0 <- ldelta$fixpar
  
  tmb_par <- list(coeff_fe_obs = matrix(unlist(Par0),nrow=ncol(X_fe_obs),ncol=1,dimnames=list(colnames(X_fe_obs))),#self$obs()$coeff_fe(), 
                  log_lambda_obs = 0, 
                  coeff_fe_hid = matrix(fixParIndex$beta0$beta,ncol=1,dimnames=list(colnames(X_fe_hid))),#self$hid()$coeff_fe(), 
                  log_lambda_hid = 0, 
                  log_delta0 = ldelta0, 
                  coeff_re_obs = 0, 
                  coeff_re_hid = 0)
  
  workBounds <- convertWorkBounds(workBounds,tmb_par,distnames)
  
  map <- NULL
  random <- NULL
  #if (is.null(S_obs)) {
    map <- c(map, list(coeff_re_obs = factor(NA), log_lambda_obs = factor(NA)))
    S_obs <- as_sparse(matrix(0, 1, 1))
    ncol_re_obs <- matrix(-1, nrow = 1, ncol = 1)
    X_re_obs <- as_sparse(rep(0, nrow(X_fe_obs)))
  #}
  #else {
  #  random <- c(random, "coeff_re_obs")
  #  tmb_par$coeff_re_obs <- rep(0, ncol(X_re_obs))
  #  tmb_par$log_lambda_obs <- log(self$lambda()$obs)
  #}
  #if (is.null(S_hid)) {
    map <- c(map, list(coeff_re_hid = factor(NA), log_lambda_hid = factor(NA)))
    S_hid <- as_sparse(matrix(0, 1, 1))
    ncol_re_hid <- matrix(-1, nrow = 1, ncol = 1)
    X_re_hid <- as_sparse(rep(0, nrow(X_fe_hid)))
  #}
  #else {
  #  random <- c(random, "coeff_re_hid")
  #  tmb_par$coeff_re_hid <- rep(0, ncol(X_re_hid))
  #  tmb_par$log_lambda_hid <- log(self$lambda()$hid)
  #}
  
  #if(is.null(betaCons) & all(is.na(fixParIndex$fixPar$beta))){
  #  fixparhid <- NULL
  #} else {
  #  if(!is.null(betaCons)){
  #    fixparhid <- betaCons
  #  } else fixparhid <- 1:length(fixParIndex$fixPar$beta)
  #  if(any(!is.na(fixParIndex$fixPar$beta))){
  #    fixparhid[which(!is.na(fixParIndex$fixPar$beta))] <- NA
  #  }
  #  names(fixparhid) <- rownames(tmb_par$coeff_fe_hid)
  #}
  #if(any(!is.na(unlist(fixParIndex$fixPar[distnames])))){
  #  fixparobs <- 1:length(unlist(fixParIndex$fixPar[distnames]))
  #  fixparobs[which(!is.na(unlist(fixParIndex$fixPar[distnames])))] <- NA
  #  names(fixparobs) <- colnames(X_fe_obs)
  #} else fixparobs <- NULL
  #fixpar <- list(hid=c(fixparhid),obs=c(fixparobs),delta0=fixParIndex$delta0)#c(self$hid()$fixpar(all = TRUE), self$obs()$fixpar(all = TRUE))
  fixpar <- list(hid=c(fixParIndex$fixPar$beta),obs=unlist(fixParIndex$fixPar[distnames]),delta0=fixParIndex$delta0)
  if(!is.null(fixpar$hid)) names(fixpar$hid) <- colnames(X_fe_hid)
  names(fixpar$obs) <- colnames(X_fe_obs)
  
  par_list <- list(coeff_fe_obs=tmb_par$coeff_fe_obs,
                   log_lambda_obs=matrix(0,0,1),
                   coeff_fe_hid=tmb_par$coeff_fe_hid,
                   log_lambda_hid=matrix(0,0,1),
                   log_delta0=tmb_par$log_delta0,
                   coeff_re_obs=matrix(0,0,1),
                   coeff_re_hid=matrix(0,0,1))#self$coeff_list()
  usernms <- c("obs", "lambda_obs", "hid", "lambda_hid", "delta0", 
               NA, NA)
  par_names <- names(par_list)
  fixpar_vec <- NULL
  for (i in seq_along(par_list)) {
    v <- par_list[[par_names[i]]]
    fix_or_not <- rep(0, length(v))
    names(fix_or_not) <- rep(par_names[i], length(v))
    fixed <- fixpar[[usernms[i]]]
    if (!is.null(fixed)) {
      tmp <- 1:length(v)
      fixed <- fixed + length(v)
      if (is.matrix(v)) {
        nms <- rownames(v)
      }
      else {
        nms <- names(v)
      }
      tmp[nms %in% names(fixed)] <- as.numeric(fixed)
      tmp <- factor(as.vector(tmp))
      ls <- list(tmp)
      names(ls) <- par_names[i]
      map <- c(map, ls)
      if (any(is.na(tmp))) {
        fix_or_not[which(is.na(tmp))] <- 1
      }
    }
    fixpar_vec <- c(fixpar_vec, fix_or_not)
  }
  coeff_array <- cbind(fixpar_vec, unlist(par_list, use.names = FALSE))
  colnames(coeff_array) <- c("fixed", "value")
  #private$coeff_array_ <- coeff_array
  statdist <- ifelse(stationary,1,0)#ifelse(self$hid()$stationary(), yes = 1, no = 0)
  
  priors <- list(coeff_fe_obs=matrix(NA,nrow(tmb_par$coeff_fe_obs),2,dimnames = list(rownames(tmb_par$coeff_fe_obs),c("mean","sd"))),
                 coeff_fe_hid=matrix(NA,nrow(tmb_par$coeff_fe_hid),2,dimnames = list(rownames(tmb_par$coeff_fe_hid),c("mean","sd"))),
                 log_delta0=matrix(NA,length(fixParIndex$delta0),2,dimnames=list(NULL,c("mean","sd"))),
                 log_lambda_obs=matrix(NA,0,2,dimnames=list(NULL,c("mean","sd"))),
                 log_lambda_hid=matrix(NA,0,2,dimnames=list(NULL,c("mean","sd"))))#self$priors()
  
  if(!is.null(prior)){
    if(!is.list(prior) || !length(names(prior))) stop("prior must be a named list")
    parCount <- 0
    for(i in distnames){
      if(!is.null(prior[[i]])){
        if(!is.matrix(prior[[i]]) || !all(dim(prior[[i]])==c(length(Par0[[i]]),2))) stop("prior$",i," must be a ",length(Par0[[i]])," x 2 matrix")
        priors$coeff_fe_obs[parCount+1:length(Par0[[i]]),] <- prior[[i]]
      }
      parCount <- parCount + length(Par0[[i]])
    }
    if(!is.null(prior$beta)){
      if(!is.matrix(prior$beta) || !all(dim(prior$beta)==c(nrow(tmb_par$coeff_fe_hid),2))) stop("prior$beta must be a ",nrow(tmb_par$coeff_fe_hid)," x 2 matrix")
      priors$coeff_fe_hid <- prior$beta
    }
  }
  
  obsvar <- obs_var(data,dist,expand=TRUE)#self$obs()$obs_var(expand = TRUE)
  datadim <- attr(obsvar, "datadim")
  
  known_states <- matrix(1, nrow = nrow(data), ncol = nbStates)
  if(!is.null(knownStates)){
    kn <- lapply(strsplit(as.character(knownStates), ","), FUN = as.numeric)
    for (i in 1:nrow(data)) {
      known_states[i, -kn[[i]]] <- 0
    }
  }
  
  tmb_dat <- list(ID = data$ID,#self$obs()$data()$ID, 
                  data = as.matrix(obsvar), 
                  datadim = datadim, known_states = known_states,#self$obs()$known_states(), 
                  n_states = nbStates, statdist = statdist, distcode = distcode, 
                  distpar = distpar, CT = ifelse(CT,1,0), delta_t = dt, dtIndex = dtIndex, kappa = kappa,
                  X_fe_obs = as_sparse(X_fe_obs), X_re_obs = as_sparse(X_re_obs), 
                  S_obs = as_sparse(S_obs), ncol_re_obs = ncol_re_obs, 
                  X_fe_hid = as_sparse(X_fe_hid), X_re_hid = as_sparse(X_re_hid), 
                  S_hid = as_sparse(S_hid), ncol_re_hid = ncol_re_hid, 
                  include_smooths = 1, ref_tpm = fixParIndex$betaRef,#self$hid()$ref(), 
                  coeff_fe_obs_prior = priors$coeff_fe_obs, 
                  coeff_fe_hid_prior = priors$coeff_fe_hid, log_lambda_obs_prior = priors$log_lambda_obs, 
                  log_lambda_hid_prior = priors$log_lambda_hid,
                  lower_fe_obs = workBounds$lower$coeff_fe_obs,
                  lower_lambda_obs = workBounds$lower$log_lambda_obs,
                  lower_fe_hid = workBounds$lower$coeff_fe_hid,
                  lower_lambda_hid = workBounds$lower$log_lambda_hid,
                  lower_delta0 = workBounds$lower$log_delta0,
                  lower_re_obs = workBounds$lower$coeff_re_obs,
                  lower_re_hid = workBounds$lower$coeff_re_hid,
                  upper_fe_obs = workBounds$upper$coeff_fe_obs,
                  upper_lambda_obs = workBounds$upper$log_lambda_obs,
                  upper_fe_hid = workBounds$upper$coeff_fe_hid,
                  upper_lambda_hid = workBounds$upper$log_lambda_hid,
                  upper_delta0 = workBounds$upper$log_delta0,
                  upper_re_obs = workBounds$upper$coeff_re_obs,
                  upper_re_hid = workBounds$upper$coeff_re_hid)
  
  #obj <- TMB::MakeADFun(c(model="momentuhmmTMB",tmb_dat), tmb_par, DLL = "momentuHMM_TMBExports", random = random, 
  #                 map = map, silent = ifelse(is.null(control$silent),FALSE,TRUE))
  #private$tmb_obj_ <- obj

  tmb_dat$include_smooths <- -1
  #private$tmb_obj_joint_ <- MakeADFun(tmb_dat, tmb_par, DLL = "hmmTMB", 
  #                                    map = map, silent = silent)
  tmb_obj <- TMB::MakeADFun(c(model="momentuhmmTMB",tmb_dat), tmb_par, DLL = "momentuHMM_TMBExports", 
                             map = map, silent = ifelse(is.null(control$silent),FALSE,control$silent))
  nllk0 <- tmb_obj$fn(tmb_obj$par)
  
  if(fit){
    
    if (is.nan(nllk0) | is.infinite(nllk0)) {
      stop(paste("Log-likelihood is NaN or infinite at starting parameters.", 
                 "Check that the data are within the domain of definition of the", 
                 "observation distributions, or try other starting parameters."))
    }
    
    control$silent <- NULL
    
    # Fit model
    args <- control
    
    if (any(c("par", "fn", "gr", "he", "hessian") %in% names(args))) {
      stop(paste("Cannot supply arguments to 'control' with name par,",
                 "fn, gr, he, or hessian. These are reserved by TMB."))
    }
    if (any(c("lower", "upper") %in% names(args))) {
      warning("'lower' and 'upper' arguments to optimx are ignored when optMethod='TMB'")
    }
    # change default method to nlminb
    tmb_obj$method <- "nlminb"
    if ("method" %in% names(args)) {
      tmb_obj$method <- args$method 
    }
    
    if (is.null(args) || !length(args)) {
      args$control <- list(kkt = FALSE, 
                           starttests = FALSE,
                           dowarn = FALSE)
    } else {
      args<- control[which(names(control) %in% c("method","itnmax","hessian"))]
      args$control <- control[which(!names(control) %in% c("hess","method","itnmax","hessian"))]
      if(is.null(args$control$kkt)) args$control$kkt <- FALSE
      if(is.null(args$control$starttests)) args$control$starttests <- FALSE
      if(is.null(args$control$dowarn)) args$control$dowarn <- FALSE
    }
    # create temporary optimization function
    opt_fn <- function(par) {as.vector(args$fn(par))}
    opt_gr <- function(par) {as.vector(args$gr(par))}
    if(isFALSE(control$hess)) {
      opt_he <- NULL
    } else {
      opt_he <- function(par) {args$he(par)}
    }
    
    # fit model 
    args <- c(tmb_obj, args)
    out <- optimx::optimx(par = args$par, 
                           fn = opt_fn, 
                           gr = opt_gr, 
                           hess = opt_he,
                           method = args$method,
                           itnmax = args$itnmax, 
                           hessian = args$hessian, 
                           control = args$control)
    best <- which.min(out$value)
    
    # Get estimates and precision matrix for all parameters
    best_par <- as.vector(out[best, 1:length(tmb_obj$par)])
    rownames(best_par) <- NULL
    names(best_par) <- names(tmb_obj$par)
    tmb_obj$par <- best_par
    tmb_obj$method <- tmb_obj$method[best]
    out <- out[best,]
    
    mod <- list()
  
    mod$minimum <- out$value
    mod$estimate <- unlist(best_par)#tmb_rep$par.fixed
    mod$code <- out$convcode
    mod$iterations <- out$fevals
    mod$tmb_obj <- tmb_obj
    if(hessian){
      tmb_rep <- TMB::sdreport(tmb_obj,
                               getJointPrecision = FALSE, 
                               skip.delta.method = FALSE)
      mod$gradient <- c(tmb_rep$gradient.fixed)
      mod$hessian <- tmb_obj$he(tmb_obj$par)
      mod$tmb_rep <- tmb_rep
    }
    if(!is.null(prior)) mod$prior <- oldpriors
    
    mod$out <- out
    
  } else {
    mod <- list()
    mod$minimum <- nllk0
    mod$estimate <- tmb_obj$par
  }
  # reorder rwdist parms
  if(!isTRUE(crwST)){
    parInd <- which(names(mod$estimate)=="coeff_fe_obs")
    convInd <- lapply(distnames,function(x) match(DMnames[[x]],names(Par0[[x]])[which(!is.na(fixParIndex$fixPar[[x]]) & !duplicated(fixParIndex$fixPar[[x]]))],nomatch=NA))
    names(convInd) <- distnames
    k <- 0
    for(i in distnames){
      if(any(is.na(convInd[[i]]))) convInd[[i]] <- convInd[[i]][-which(is.na(convInd[[i]]))]
      if(!all(is.na(fixParIndex$fixPar[[i]]))){
        maxf <- max(convInd[[i]],na.rm=TRUE)
        convInd[[i]] <- k + convInd[[i]]
        k <- k + maxf
      }
    }
    convInd <- unlist(convInd)
    mod$estimate[parInd] <- mod$estimate[convInd]
    if(hessian){
      if(max(convInd)<length(mod$estimate)){
        mod$hessian <- mod$hessian[c(convInd,(max(convInd)+1):length(mod$estimate)),c(convInd,(max(convInd)+1):length(mod$estimate))]
        mod$gradient <- mod$gradient[c(convInd,(max(convInd)+1):length(mod$estimate))]
      } else {
        mod$hessian <- mod$hessian[convInd,convInd]
        mod$gradient <- mod$gradient[convInd]
      }
    }
  }
  
  return(mod)
}

dist_code <- function(dist,zeroInflation,oneInflation) {
  
  if(grepl("cat",dist)){
    dimCat <- as.integer(gsub("cat","",dist))
    if(is.na(dimCat)) stop("categorical distributions must be specified using paste0('cat',k), where k is the number of categories (e.g. 'cat3', 'cat12', etc.)")
    if(dimCat<2) stop("categorical distribution must have at least 2 categories")
    dist <- "cat"
  } else if(grepl("ctds",dist)){
    dimCat <- attr(dist,"directions") + 1
  }
  
  switch(dist,
         "bern"={code <- 1},
         "beta"={ 
           if(zeroInflation | oneInflation) code <- NULL
           else code <- 0},
         "cat"={code <- 2},
         "crwrice"={
           if(zeroInflation) code <- NULL
           else code <- 29},
         "crwvm"={code <- 30},
         "ctcrw"={code <- 31},
         "ctds"={code <- NULL},
         "exp"={
           if(zeroInflation) code <- NULL
           else code <- 4},
         "gamma"={
           if(zeroInflation) code <- 21
           else code <-  7},
         "logis"={code <- NULL},
         "lnorm"={if(zeroInflation) code <- NULL
         else code <- 8},
         "negbinom"={code <- 10},
         "norm"={code <- 11},
         "mvnorm2"={code <- 9},
         "mvnorm3"={code <- 9},
         "rw_norm"={code <- 11},
         "rw_mvnorm2"={code <- 9},
         "rw_mvnorm3"={code <- 9},
         "pois"={code <- 12},
         "t"={code <- 13},
         "vm"={code <- 16},
         "vmConsensus"={code <- NULL},
         "weibull"={
           if(zeroInflation) code <- NULL
           else code <- 17},
         "wrpcauchy"={code <- 18}
  )
  if(is.null(code)) stop("Sorry,",ifelse(zeroInflation," zero-inflated ",ifelse(oneInflation," one-inflated ","")) ,"'",dist,"' data stream distributions are not currently supported with TMB ")
  else return(code)
}

find_re <- function(form) {
  term_labs <- attr(stats::terms(form), "term.labels")
  pattern <- "^s\\((.*), bs = \"re\"\\)$"
  var_names <- gsub(pattern = pattern, "\\1", term_labs)
  which_re <- grep(pattern = pattern, term_labs)
  var_re <- var_names[which_re]
  if (length(var_re) > 0) {
    var_re <- strsplit(var_re, split = ",")[[1]][1]
  }
  return(var_re)
}

make_matrices <- function (formulas, data, new_data = NULL) 
{
  X_list_fe <- list()
  X_list_re <- list()
  S_list <- list()
  ncol_fe <- NULL
  ncol_re <- NULL
  names_fe <- NULL
  names_re <- NULL
  names_ncol_re <- NULL
  start <- 1
  forms <- unlist(formulas)
  names <- names(forms)
  for (k in seq_along(forms)) {
    form <- forms[[k]]
    var_re <- find_re(form)
    for (var in var_re) {
      if (!inherits(data[[var]], "factor")) {
        data[[var]] <- factor(data[[var]])
        warning(paste0("'", var, "' is included as a random effect but is ", 
                       "not a factor - changing to factor."))
      }
    }
    if (is.null(new_data)) {
      gam_setup <- mgcv::gam(formula = stats::update(form, dummy ~ .), 
                       data = cbind(dummy = 1, data), fit = FALSE)
      Xmat <- gam_setup$X
      term_names <- gam_setup$term.names
    }
    else {
      gam_setup0 <- mgcv::gam(formula = stats::update(form, dummy ~ 
                                           .), data = cbind(dummy = 1, data))
      gam_setup <- mgcv::gam(formula = stats::update(form, dummy ~ .), 
                       data = cbind(dummy = 1, data), fit = FALSE)
      Xmat <- stats::predict(gam_setup0, newdata = new_data, type = "lpmatrix")
      term_names <- gam_setup$term.names
    }
    X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
    subnames_fe <- paste0(names[k], ".", term_names[1:gam_setup$nsdf])
    names_fe <- c(names_fe, subnames_fe)
    X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
    if (ncol(X_list_re[[k]]) > 0) {
      subnames_re <- paste0(names[k], ".", term_names[-(1:gam_setup$nsdf)])
      names_re <- c(names_re, subnames_re)
    }
    S_list[[k]] <- bdiag_check(gam_setup$S)
    ncol_fe <- c(ncol_fe, gam_setup$nsdf)
    if (length(gam_setup$smooth) > 0) {
      sub_ncol_re <- matrix(1, nrow = 2, ncol = length(gam_setup$S))
      colnames(sub_ncol_re) <- 1:ncol(sub_ncol_re)
      start_s <- 1
      for (s in 1:length(gam_setup$smooth)) {
        npen <- length(gam_setup$smooth[[s]]$S)
        npar <- ncol(gam_setup$smooth[[s]]$S[[1]])
        sub_ncol_re[, (start_s:(start_s + npen - 1))] <- c(start, 
                                                           start + npar - 1)
        colnames(sub_ncol_re)[start_s:(start_s + npen - 
                                         1)] <- rep(gam_setup$smooth[[s]]$label, npen)
        s_terms <- gsub("(.*)\\..*", "\\1", names_re[sub_ncol_re[1, 
                                                                 s]:sub_ncol_re[2, s]])
        names_ncol_re <- c(names_ncol_re, rep(unique(s_terms), 
                                              npen))
        start <- start + npar
        start_s <- start_s + npen
      }
      ncol_re <- cbind(ncol_re, sub_ncol_re)
    }
  }
  colnames(ncol_re) <- names_ncol_re
  X_fe <- bdiag_check(X_list_fe)
  colnames(X_fe) <- names_fe
  X_re <- bdiag_check(X_list_re)
  colnames(X_re) <- names_re
  S <- bdiag_check(S_list)
  return(list(X_fe = X_fe, X_re = X_re, S = S, X_list_fe = X_list_fe, 
              X_list_re = X_list_re, S_list = S_list, ncol_fe = ncol_fe, 
              ncol_re = ncol_re))
}
as_sparse <- function(x) {
  if (length(dim(x)) < 2) {
    x <- matrix(x, ncol = 1)
  }
  mat <- suppressMessages(methods::as(x, "dgTMatrix"))
  return(mat)
}

obs_var <- function(data,dist,expand = FALSE) {
  distnames <- names(dist)
  if("step" %in% distnames) steps <- data$step
  if("angle" %in% distnames) angles <- data$angle
  for(i in distnames){
    if(dist[[i]]=="mvnorm2" || dist[[i]]=="rw_mvnorm2"){
      data[[i]] <- split(as.matrix(data[,c(paste0(i,".x"),paste0(i,".y"))]),seq(nrow(data)))
      data[[paste0(i,".x")]] <- data[[paste0(i,".y")]] <- NULL
    } else if(dist[[i]]=="mvnorm3" || dist[[i]]=="rw_mvnorm3"){
      data[[i]] <- split(as.matrix(data[,c(paste0(i,".x"),paste0(i,".y"),paste0(i,".z"))]),seq(nrow(data)))
      data[[paste0(i,".x")]] <- data[[paste0(i,".y")]]<- data[[paste0(i,".z")]] <- NULL
    } else if(dist[[i]]=="crwrice" | dist[[i]]=="crwvm"){
      if(dist[[i]]=="crwrice"){
        genData <- steps
        step_tm1 <- c(0,steps[-length(steps)])
        step_tm1[c(1,1+cumsum(table(data$ID)[-(length(unique(data$ID)))]))] <- .Machine$double.xmin
        genInd <- which(!is.na(genData))
        step_tm1[-genInd] <- NA
        data[[i]] <- split(cbind(step_tm1,genData),seq(nrow(data)))
      } else if(dist[[i]]=="crwvm"){
        genData <- angles
        step_tm1 <- c(0,steps[-length(steps)])
        step_tm1[c(1,1+cumsum(table(data$ID)[-(length(unique(data$ID)))]))] <- .Machine$double.xmin
        step <- steps
        genInd <- which(!is.na(genData))
        step_tm1[-genInd] <- NA
        step[-genInd] <- NA
        data[[i]] <- split(cbind(step_tm1,step,genData),seq(nrow(data)))
      }
    } else if(dist[[i]]=="ctcrw"){
      tmpdat <- data[,c(paste0(i,".x"),paste0(i,".y"),paste0(i,".x_tm1"),paste0(i,".y_tm1"))]
      tmpdat[[paste0(i,".x_tm2")]] <- 0
      tmpdat[[paste0(i,".y_tm2")]] <- 0
      tmpdat$dtm1 <- 1
      for(zoo in unique(data$ID)){
        tInd <- which(data$ID==zoo)
        tmpdat[[paste0(i,".x_tm2")]][tInd] <- dplyr::lag(tmpdat[[paste0(i,".x_tm1")]][tInd],n=1,default = 0)
        tmpdat[[paste0(i,".y_tm2")]][tInd] <- dplyr::lag(tmpdat[[paste0(i,".y_tm1")]][tInd],n=1,default = 0)
        tmpdat[[paste0(i,".x_tm2")]][tInd][1] <- tmpdat[[paste0(i,".x_tm1")]][tInd[1]] # no crw term for first observation
        tmpdat[[paste0(i,".y_tm2")]][tInd][1] <- tmpdat[[paste0(i,".y_tm1")]][tInd[1]] # no crw term for first observation
        if(!is.null(data$dt)) tmpdat$dtm1[tInd] <- dplyr::lag(data$dt[tInd],n=1,default = 1)
      }
      data[[i]] <- split(as.matrix(tmpdat),seq(nrow(data)))
      data[[paste0(i,".x")]] <- data[[paste0(i,".y")]] <- NULL
      data[[paste0(i,".x_tm1")]] <- data[[paste0(i,".y_tm1")]] <- NULL
    }
  }
  obs_names <- distnames
  obs_var <- data[, obs_names, drop = FALSE]
  datadim <- rep(1, ncol(obs_var))
  if (expand) {
    multivar <- sapply(obs_var, is.list)
    if (any(multivar)) {
      wh <- which(multivar)
      for (i in 1:length(wh)) {
        if(i==1) ind <- wh[i]
        else ind <- sum(datadim[1:(wh[i]-1)])+1
        v <- do.call(rbind, obs_var[[ind]])
        datadim[wh[i]] <- ncol(v)
        tmp <- NULL
        tmpnms <- NULL
        if (wh[i] > 1) {
          tmp <- obs_var[,1:(ind - 1)]
          tmpnms <- c(tmpnms, colnames(obs_var)[1:(ind - 1)])
        }
        tmp <- cbind(tmp, v)
        tmpnms <- c(tmpnms, rep(names(wh[i]), ncol(v)))
        if (ind < ncol(obs_var)) {
          tmp <- cbind(tmp, obs_var[,(ind + 1):ncol(obs_var),drop=FALSE])
          tmpnms <- c(tmpnms, colnames(obs_var)[(ind + 1):ncol(obs_var)]) 
        }
        obs_var <- tmp
        colnames(obs_var) <- tmpnms
      }
    }
  }
  attributes(obs_var)$datadim <- datadim
  return(obs_var)
}

bdiag_check <- function(...) {
  args <- list(...)[[1]]
  check <- sapply(args, function(arg) {
    !is.null(dim(arg)) | length(arg) > 0
  })
  if (length(check) == 0) 
    return(NULL)
  else return(Matrix::bdiag(args[check]))
}

#expandDM <- function(nbObs,DMinputs){
#  DM <- DMinputs$fullDM
#  npar <- sum(unlist(lapply(DM,ncol)))
#  tmbDM <- matrix(0,nrow=nbObs*npar,ncol=npar,dimnames=list(NULL,unlist(lapply(DM,colnames))))
#  parCount <- 1
#  for(i in names(DM)){
#    if(DMinputs$DMind[[i]]){
#      for(j in 1:ncol(DM[[i]])){
#        tmbDM[(parCount-1)*nbObs+1:nbObs,parCount] <- rep(DM[[i]][j,j],nbObs)
#        parCount <- parCount + 1
#      }
#    } else {
#      
#    }
#  }
#  return(Matrix::bdiag(tmbDM))
#}

make_formulas <- function(input_forms, var_names, par_names, nbStates){
  output_forms <- list()
  form_names <- names(input_forms)
  for (i in 1:length(var_names)) {
    mch <- match(var_names[i], form_names)
    if (!is.na(mch)) {
      var_forms <- input_forms[[mch]]
    }
    else {
      var_forms <- NULL
    }
    var_forms_new <- list()
    for (j in 1:length(par_names[[i]])) {
      par_mch <- match(par_names[[i]][[j]], names(var_forms))
      if (!is.na(par_mch)) {
        form <- var_forms[[par_mch]]
      }
      else {
        form <- as.formula("~ 1")
      }
      form_terms <- terms(form, specials = paste0("state", 
                                                  1:nbStates),keep.order=TRUE)
      labs <- attr(form_terms, "term.labels")
      covs <- gsub(pattern = "^state[0-9]*\\((.*)\\)$", 
                   replacement = "\\1", x = labs)
      which_all_states <- which(!seq_along(labs) %in% 
                                  match(rownames(attr(form_terms,"factors"))[unlist(attr(form_terms,"specials"))],colnames(attr(form_terms,"factors"))))
      state_forms <- list()
      for (s in 1:nbStates) {
        which_this_state <- match(rownames(attr(form_terms,"factors"))[attr(form_terms,"specials")[[paste0("state",s)]]],colnames(attr(form_terms,"factors")))
        if (attr(form_terms, "intercept") == 1) {
          new_form <- "~ 1"
        }
        else {
          new_form <- "~ -1"
        }
        for (k in which_all_states) new_form <- paste0(new_form, 
                                                       " + ", covs[k])
        for (k in which_this_state) new_form <- paste0(new_form, 
                                                       " + ", covs[k])
        state_forms[[paste0("state", s)]] <- stateFormulas(as.formula(new_form),nbStates=1)[[1]]
      }
      var_forms_new[[j]] <- state_forms
      names(var_forms_new)[j] <- par_names[[i]][[j]]
    }
    output_forms[[i]] <- var_forms_new
    names(output_forms)[i] <- var_names[i]
  }
  return(output_forms)
}

getTPMmat <- function(formula,nbStates,betaRef,data){
  
  formula <- matrix(deparse1(formula), 
                    nrow = nbStates, 
                    ncol = nbStates)
  
  ref_mat <- matrix(0,nbStates,nbStates)
  for(j in 1:nbStates){
    ref_mat[j,betaRef[j]] <- 1
  }
  
  ls_form_char <- as.list(t(formula)[!t(ref_mat)])
  ls_form <- lapply(ls_form_char, function(form_char) {
    if(form_char == ".")
      return(as.formula("~1"))
    else
      return(as.formula(form_char))
  })
  
  # Names for transition probabilities
  tr_names <- paste0("S", rep(1:nbStates, each = nbStates), 
                     ">S", rep(1:nbStates, nbStates))
  names(ls_form) <- tr_names[-which(t(ref_mat) == 1)]
  return(make_matrices(ls_form,data))
}

get_tmbDelta <- function(delta0,nbStates,nbAnimals,stationary,formulaDelta,covsDelta,fixPar,workBounds){
  ldelta0 <- matrix(0,nbAnimals*(nbStates-1),1,dimnames=list(paste0("ID:",rep(1:nbAnimals,nbStates-1),".state",rep(2:(nbStates),each=nbAnimals))))#self$hid()$delta0(log = TRUE, as_matrix = FALSE)
  fixpar <- NULL
  if(!stationary){
    if(is.null(formulaDelta)){
      ldelta0 <- matrix(rep(delta0,each=nbAnimals),nbAnimals*(nbStates-1),1,dimnames=list(paste0("ID:",rep(1:nbAnimals,nbStates-1),".state",rep(2:(nbStates),each=nbAnimals))))
      fixPar <- rep(fixPar,each=nbAnimals)
      workBounds <- matrix(rep(workBounds,each=nbAnimals),nrow=nbAnimals*(nbStates-1),ncol=2)
      names(fixPar) <- rownames(ldelta0)
    } else {
      dterms <- stats::terms(formulaDelta)
      if(length(attr(dterms,"factors"))){
        if(attr(dterms,"intercept")) stop("formulaDelta cannot include intercept term when fitting using TMB. For example, use '~0+ID' instead of '~ID'")
        if("ID" %in% attr(dterms,"term.labels")){
          ldelta0 <- matrix(delta0,nbAnimals*(nbStates-1),1,dimnames=list(paste0("ID:",rep(1:nbAnimals,nbStates-1),".state",rep(2:(nbStates),each=nbAnimals))))
          names(fixPar) <- rownames(ldelta0)
        } else {
          ## need to possibly add other factor covariates by adjusting ldelta0, fixpar, and map (and adjust momentuHMM delta parameters) accordingly based on unique combinations for each ID
          stop("formulaDelta = '",formulaDelta,"'  is not currently supported in models fitted by TMB. formulaDelta must be 'NULL', '~1', or '~0+ID'")
        }
      } else {
        ldelta0 <- matrix(delta0,nbAnimals*(nbStates-1),1,dimnames=list(paste0("ID:",rep(1:nbAnimals,nbStates-1),".state",rep(2:(nbStates),each=nbAnimals))))
        fixPar <- rep(fixPar,each=nbAnimals)
        workBounds <- matrix(workBounds,nrow=nbAnimals*(nbStates-1),ncol=2,byrow=TRUE)
        names(fixPar) <- rownames(ldelta0)
      }
    }
  } else {
    fixPar <- rep.int(NA,length(ldelta0))
    names(fixPar) <- rownames(ldelta0)
    workBounds <- matrix(c(-Inf,Inf),nbAnimals*(nbStates-1),2,byrow=TRUE)
  }
  return(list(ldelta0=ldelta0,fixpar=fixPar,workBounds=workBounds))
}

convertWorkBounds <- function(workBounds,par,distnames){ 
  
  parnames <- names(par)
  lower <- lapply(par,function(x) rep(-Inf,length(x)))
  upper <- lapply(par,function(x) rep(Inf,length(x)))
  
  lower$coeff_fe_obs <- unlist(lapply(workBounds[distnames],function(x) x[,1]))
  upper$coeff_fe_obs <- unlist(lapply(workBounds[distnames],function(x) x[,2]))
  
  lower$coeff_fe_hid <- workBounds$beta[,1]
  upper$coeff_fe_hid <- workBounds$beta[,2]
  
  lower$log_delta0 <- workBounds$delta[,1]
  upper$log_delta0 <- workBounds$delta[,2]
  
  return(list(lower=lower,upper=upper))
}

reorderFix <- function(fixPar){
  tmp <- c(fixPar)
  if(!all(is.na(tmp))){
    if(!all(tmp[which(!duplicated(tmp) & !is.na(tmp))]==(1:length(tmp[which(!duplicated(tmp) & !is.na(tmp))])))){
      k <- 1
      newtmp <- tmp
      for(j in unique(tmp)){
        if(!is.na(j)){
          newtmp[which(tmp==j)] <- k
          k <- k + 1
        }
      }
      fixPar <- newtmp
    } 
  }
  return(fixPar)
}
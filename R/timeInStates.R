
#' Calculate proportion of time steps assigned to each state (i.e. \dQuote{activity budgets})
#' 
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{HMMfits}} object.
#' @param by A character vector indicating any groupings by which to calculate the proportions, such as individual (\dQuote{ID}) or group-level (e.g. sex or age class) covariates. Default is \code{NULL} (no groupings are used).
#' @param alpha Significance level for calculating confidence intervals of pooled estimates. Default: 0.95. Ignored unless \code{m} is a \code{\link{miHMM}} or \code{\link{HMMfits}} object. 
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing). Ignored unless \code{m} is a \code{\link{miHMM}} or \code{\link{HMMfits}} object. 
#' 
#' @return If \code{m} is a \code{\link{momentuHMM}} object, a data frame containing the estimated activity budgets for each state (grouped according to \code{by}).  If \code{m} is a \code{\link{miHMM}} or \code{\link{HMMfits}} object, a list containing the activity budget
#' estimates, standard errors, lower bounds, and upper bounds across all imputations.
#' 
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' timeInStates(m)
#' timeInStates(m, by = "ID")
#' 
#' @export
timeInStates <- function(m, by = NULL, alpha = 0.95, ncores = 1) {
  UseMethod("timeInStates")
}

#' @method timeInStates momentuHMM
#' @export
#' @rdname timeInStates
timeInStates.momentuHMM <- function(m, by = NULL, alpha = 0.95, ncores = 1){
  
  if(any(!(by %in% names(m$data)))) stop(by[which(!(by %in% names(m$data)))]," not found in data")
  
  nbStates <- length(m$stateNames)
  if(nbStates>1) states <- momentuHMM::viterbi(m)
  else states <- rep(1,nrow(m$data))
  
  get_combins(m,states,by)
}

#' @method timeInStates HMMfits
#' @export
#' @rdname timeInStates
#' @importFrom magrittr %>%
#' @importFrom dplyr funs group_by_at summarise_at
timeInStates.HMMfits <- function(m, by = NULL, alpha = 0.95, ncores = 1){
  
  i <- . <- NULL #gets rid of no visible binding for global variable 'i' and '.' NOTE in R cmd check
  
  simind <- unlist(lapply(m,is.momentuHMM))
  nsims <- length(which(simind))
  if(nsims<1) stop("'HMMfits' must be a list comprised of momentuHMM objects")
  checkmove <- which(!simind)
  if(length(checkmove)) {
    m[checkmove]<-NULL
    warning("The following imputations are not momentuHMM objects and will be ignored: ",paste(checkmove,collapse=", "))
  }
  nsims <- length(m)
  
  stateNames <- m[[1]]$stateNames
  nbStates <- length(stateNames)
  if(any(!(by %in% names(m[[1]]$data)))) stop(by[which(!(by %in% names(m[[1]]$data)))]," not found in data")
  
  cat("Decoding state sequences and probabilities for each imputation... ")
  registerDoParallel(cores=ncores)
  im_states <- foreach(i = 1:nsims) %dopar% {momentuHMM::viterbi(m[[i]])}
  stopImplicitCluster()
  cat("DONE\n")
  
  registerDoParallel(cores=ncores)
  xmat <- foreach(i = 1:nsims, .combine = rbind) %dopar% {
    get_combins(m[[i]],im_states[[i]],by)
  }
  stopImplicitCluster()
  
  if(!is.null(by)){
    
    n <- as.data.frame(
      xmat %>%
        group_by_at(by) %>%
        summarise_at(stateNames, funs(sum(!is.na(.)))))
    
    if(any(n[,stateNames]<2)) warning("need at least 2 simulations for each 'by' combination with valid point and variance estimates")
    
    xbar <- as.data.frame(
      xmat %>%
          group_by_at(by) %>%
          summarise_at(stateNames, funs(mean(., na.rm=TRUE))))
  
    MI_se <- as.data.frame(
      xmat %>%
        group_by_at(by) %>%
        summarise_at(stateNames, funs(sqrt((sum(!is.na(.))+1)/sum(!is.na(.)) * var(., na.rm=TRUE)))))
    
    lower <- as.data.frame(
      xmat %>%
        group_by_at(by) %>%
        summarise_at(stateNames, funs(probCI(mean(., na.rm=TRUE),sqrt((sum(!is.na(.))+1)/sum(!is.na(.)) * var(., na.rm=TRUE)),qt(1-(1-alpha)/2,df=sum(!is.na(.))-1),"lower"))))
    
    upper <- as.data.frame(
      xmat %>%
        group_by_at(by) %>%
        summarise_at(stateNames, funs(probCI(mean(., na.rm=TRUE),sqrt((sum(!is.na(.))+1)/sum(!is.na(.)) * var(., na.rm=TRUE)),qt(1-(1-alpha)/2,df=sum(!is.na(.))-1),"upper"))))
    
    combins <- list(est=xbar,se=MI_se,lower=lower,upper=upper)
  
  } else {
    
    n <- apply(!is.na(xmat),2,sum)
    
    if(any(n<2)) warning("need at least 2 simulations with valid point and variance estimates")
    
    xbar <- apply(xmat,2,mean,na.rm=TRUE)
    B_m <- apply(xmat,2,var,na.rm=TRUE)
    
    MI_se <- sqrt((n+1)/n * B_m)
    
    dfs <- n-1
    quantSup <- qt(1-(1-alpha)/2,df=dfs)
    
    lower <- probCI(xbar,MI_se,quantSup,"lower")
    upper <- probCI(xbar,MI_se,quantSup,"upper")
    
    combins <- list(est=xbar,se=MI_se,lower=lower,upper=upper)
    names(combins$est) <- stateNames
    names(combins$se) <- stateNames
    names(combins$lower) <- stateNames
    names(combins$upper) <- stateNames
  }
  
  combins
}

#' @method timeInStates miHMM
#' @export
#' @rdname timeInStates
timeInStates.miHMM <- function(m, by = NULL, alpha = 0.95, ncores = 1){
  timeInStates(m$HMMfits, by = by, alpha = alpha, ncores = ncores)
}

get_combins <- function(m,states,by){
  
  nbStates <- length(m$stateNames)
  
  if(any(!(by %in% names(m$data)))) stop(by[which(!(by %in% names(m$data)))]," not found in data")
  
  if(length(by)){
    combins <- base::expand.grid(lapply(m$data[by],unique))
    combins <- cbind(combins,matrix(0,nrow(combins),nbStates))
    for(i in 1:nrow(combins)){
      counts<-hist(states[apply(m$data[by]==matrix(rep(unlist(combins[i,by]),each=nrow(m$data))),1,all)],breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts
      time <- counts/sum(counts)
      combins[i,ncol(combins)-(nbStates-1):0] <- time
    }
    colnames(combins)[ncol(combins)-(nbStates-1):0] <- m$stateNames
  } else {
    counts<-hist(states,breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts
    combins <- counts/sum(counts)
    combins <- as.data.frame(matrix(combins,1,dimnames=list(NULL,m$stateNames)))
  }
  combins
}


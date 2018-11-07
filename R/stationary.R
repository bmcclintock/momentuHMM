
#' Stationary state probabilities
#'
#' Calculates the stationary probabilities of each state based on
#' covariate values.
#'
#' @param model \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object
#' @param covs Either a data frame or a design matrix of covariates. If \code{covs} is not provided, then the stationary probabilties are calculated based on the covariate data for each time step.
#'
#' @return A list of length \code{model$conditions$mixtures} where each element is a matrix of stationary state probabilities for each mixture. For each matrix, each row corresponds to
#' a row of covs, and each column corresponds to a state.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' # data frame of covariates
#' stationary(m, covs = data.frame(cov1 = 0, cov2 = 0))
#'
#' # design matrix (each column corresponds to row of m$mle$beta)
#' stationary(m, covs = matrix(c(1,0,cos(0)),1,3))
#'
#' @export
#'
stationary <- function (model, covs) {
  UseMethod("stationary")
}

#' @method stationary momentuHMM
#' @export
stationary.momentuHMM <- function(model, covs)
{
    model <- delta_bc(model)
    
    nbStates <- length(model$stateNames)
    beta <- model$mle$beta

    if(nbStates==1)
        stop("No state probabilities (1-state model).")

    formula<-model$conditions$formula
    newForm <- newFormulas(formula,nbStates)
    formulaStates <- newForm$formulaStates
    formterms <- newForm$formterms
    newformula <- newForm$newformula
    recharge <- newForm$recharge
    
    if(missing(covs)){
      covs <- model$rawCovs
      if(length(covs)){
        tmpSplineInputs<-getSplineFormula(newformula,model$data,covs)
      } else {
        tmpSplineInputs<-getSplineFormula(newformula,model$data,model$data[1,])
      }
      covMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
    } else if(is.data.frame(covs)){
      if(is.null(recharge))
        if(!all(names(covs) %in% names(model$data))) stop('invalid covs specified')
      else 
        if(!all(names(covs) %in% c(names(model$data),"recharge"))) stop('invalid covs specified')
      if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(model$data$ID))
      for(j in names(model$data)[which(names(model$data) %in% names(covs))]){
        if(inherits(model$data[[j]],"factor")) covs[[j]] <- factor(covs[[j]],levels=levels(model$data[[j]]))
        if(any(is.na(covs[[j]]))) stop("check value(s) for ",j)
      }
      
      tmpSplineInputs<-getSplineFormula(newformula,model$data,covs)
      
      covMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
    } else if(is.matrix(covs)){
      covMat <- covs
    } else stop("covs must either be a data frame or a matrix")
    
    aInd <- NULL
    nbAnimals <- length(unique(model$data$ID))
    for(i in 1:nbAnimals){
      aInd <- c(aInd,which(model$data$ID==unique(model$data$ID)[i])[1])
    }
    
    if(!is.null(recharge)){
      if(is.matrix(covs)) stop("covs must be provided as a data frame for recharge models")
      g0covs <- model.matrix(recharge$g0,model$data[aInd,])
      nbG0covs <- ncol(g0covs)-1
      recovs <- model.matrix(recharge$theta,model$data)
      nbRecovs <- ncol(recovs)-1
      model$data$recharge<-rep(0,nrow(model$data))
      for(i in 1:nbAnimals){
        idInd <- which(model$data$ID==unique(model$data$ID)[i])
        if(nbRecovs){
          g0 <- model$mle$g0 %*% t(g0covs[i,,drop=FALSE])
          theta <- model$mle$theta
          model$data$recharge[idInd] <- cumsum(c(g0,theta%*%t(recovs[idInd[-length(idInd)],])))
        }
      }
      if(!length(covs)) covs <- model$data
      newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
      tmpSplineInputs<-getSplineFormula(newformula,model$data,covs)
      # check that all covariates are provided
      ck1 <- tryCatch(model.matrix(recharge$theta,tmpSplineInputs$covs),error=function(e) e)
      ck2 <- tryCatch(model.matrix(tmpSplineInputs$formula,tmpSplineInputs$covs),error=function(e) e)
      if(inherits(ck1,"error")) stop("covs not specified correctly -- ",ck1)
      if(inherits(ck2,"error")) stop("covs not specified correctly -- ",ck2)
    }

    mixtures <- model$conditions$mixtures
  
    probs <- list()
    
    for(mix in 1:mixtures){
      # all transition matrices
      tmpbeta <- beta[(mix-1)*ncol(covMat)+1:ncol(covMat),,drop=FALSE]
      if(is.null(recharge))
        allMat <- trMatrix_rcpp(nbStates=nbStates, beta=tmpbeta, covs=covMat, betaRef=model$conditions$betaRef)
      else {
        gamInd<-(length(model$mod$estimate)-(nrow(tmpbeta))*nbStates*(nbStates-1)*mixtures+1):(length(model$mod$estimate))-(mixtures-1)-ifelse(nbRecovs,(nbRecovs+1)+(nbG0covs+1),0)-ncol(model$covsDelta)*(nbStates-1)*(!model$conditions$stationary)*mixtures
        allMat <- array(unlist(lapply(split(tmpSplineInputs$covs,1:nrow(covs)),function(x) tryCatch(get_gamma_recharge(model$mod$estimate[c(gamInd[unique(c(model$conditions$betaCons))],length(model$mod$estimate)-nbRecovs:0)],covs=x,formula=tmpSplineInputs$formula,recharge=recharge,nbStates=nbStates,betaRef=model$conditions$betaRef,betaCons=model$conditions$betaCons,workBounds=rbind(model$conditions$workBounds$beta,model$conditions$workBounds$theta),mixture=mix),error=function(e) NA))),dim=c(nbStates,nbStates,nrow(covs)))
      }
      
      tryCatch({
          # for each transition matrix, derive corresponding stationary distribution
          probs[[mix]] <- apply(allMat, 3,
                         function(gamma)
                             solve(t(diag(nbStates)-gamma+1),rep(1,nbStates)))
          probs[[mix]] <- t(probs[[mix]])
      },
      error = function(e) {
          stop(paste("The stationary probabilities cannot be calculated",
                     "for these covariate values (singular system)."))
      })
  
      colnames(probs[[mix]]) <- model$stateNames
    }
    return(probs)
}

#' @method stationary miSum
#' @export
stationary.miSum <- function(model, covs)
{
  model$mle <- lapply(model$Par$real,function(x) x$est)
  model$mle$beta <- model$Par$beta$beta$est
  model$mod <- list()

  stationary(momentuHMM(model),covs)
}

#' @method stationary miHMM
#' @export
stationary.miHMM <- function(model, covs)
{
  stationary(model$miSum,covs)
}

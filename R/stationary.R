
#' Stationary state probabilities
#'
#' Calculates the stationary probabilities of each state based on
#' covariate values.
#'
#' @param model \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object
#' @param covs Either a data frame or a design matrix of covariates. If \code{covs} is not provided, then the stationary probabilties are calculated based on the covariate data for each time step.
#'
#' @return Matrix of stationary state probabilities. Each row corresponds to
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

    if(missing(covs)){
        covs <- model$rawCovs
        if(length(covs)){
            tmpSplineInputs<-getSplineFormula(newformula,model$data,covs)
        } else {
            tmpSplineInputs<-getSplineFormula(newformula,model$data,model$data[1,])
        }
        covMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
    } else if(is.data.frame(covs)){
        if(!all(names(covs) %in% names(model$data))) stop('invalid covs specified')
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

    # all transition matrices
    allMat <- trMatrix_rcpp(nbStates=nbStates, beta=beta, covs=covMat, betaRef=model$conditions$betaRef)

    tryCatch({
        # for each transition matrix, derive corresponding stationary distribution
        probs <- apply(allMat, 3,
                       function(gamma)
                           solve(t(diag(nbStates)-gamma+1),rep(1,nbStates)))
        probs <- t(probs)
    },
    error = function(e) {
        stop(paste("The stationary probabilities cannot be calculated",
                   "for these covariate values (singular system)."))
    })

    colnames(probs) <- model$stateNames

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

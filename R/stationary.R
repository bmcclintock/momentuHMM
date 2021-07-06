
#' Stationary state probabilities
#'
#' Calculates the stationary probabilities of each state based on
#' covariate values.
#'
#' @param model \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object
#' @param covs Either a data frame or a design matrix of covariates. If \code{covs} is not provided, then the stationary probabilties are calculated based on the covariate data for each time step.
#' @param covIndex Integer vector indicating specific rows of the data to be used in the calculations. This can be useful for reducing unnecessarily long computation times, e.g., when \code{formula} includes factor covariates (such as \code{ID}) but no temporal covariates. Ignored unless \code{covs} is missing.
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
#' # get stationary distribution for first 3 observations
#' stationary(m, covIndex = c(1,2,3))
#'
#' @export
#'
stationary <- function (model, covs, covIndex) {
  UseMethod("stationary")
}

#' @method stationary momentuHMM
#' @export
stationary.momentuHMM <- function(model, covs, covIndex = NULL)
{
    model <- delta_bc(model)
    
    nbStates <- length(model$stateNames)
    beta <- model$mle$beta

    if(nbStates==1)
        stop("No state probabilities (1-state model).")
    
    if(inherits(model,"hierarchical")) installDataTree()

    formula<-model$conditions$formula
    newForm <- newFormulas(formula,nbStates,model$conditions$betaRef,hierarchical=TRUE)
    newformula <- newForm$newformula
    recharge <- newForm$recharge
    
    if(!missing(covs)){
      if(!is.null(covIndex)) stop("Either 'covs' or 'covIndex' can be specified (not both)")
      if(is.data.frame(covs)){
        if(is.null(recharge)){
          if(!all(names(covs) %in% names(model$data))) stop('invalid covs specified')
        } else {
          if(!inherits(model,"hierarchical")) if(!all(names(covs) %in% c(names(model$data),"recharge"))) stop('invalid covs specified')
          else if(!all(names(covs) %in% c(names(model$data),paste0("recharge",levels(model$data$level))))) stop('invalid covs specified')
        }
        if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(model$data$ID))
        if(inherits(model,"hierarchical") && any(names(covs) %in% "level")) stop("covs$level cannot be specified for hierarchical models")
        for(j in names(model$data)[which(names(model$data) %in% names(covs))]){
          if(inherits(model$data[[j]],"factor")) covs[[j]] <- factor(covs[[j]],levels=levels(model$data[[j]]))
          if(any(is.na(covs[[j]]))) stop("check value(s) for ",j)
        }
        
        tmpSplineInputs<-getSplineFormula(newformula,model$data,covs)
  
      } else if(is.matrix(covs)){
        covMat <- covs
      } else stop("covs must either be a data frame or a matrix")
    }
    
    aInd <- NULL
    nbAnimals <- length(unique(model$data$ID))
    for(i in 1:nbAnimals){
      aInd <- c(aInd,which(model$data$ID==unique(model$data$ID)[i])[1])
    }
    
    if(!is.null(recharge)){
      reForm <- formatRecharge(nbStates,formula,model$conditions$betaRef,data=model$data,par=list(g0=model$mle$g0,theta=model$mle$theta))
      newformula <- reForm$newformula
      recharge <- reForm$recharge
      hierRecharge <- reForm$hierRecharge
      model$data[colnames(reForm$newdata)] <- reForm$newdata
      g0covs <- reForm$g0covs
      nbG0covs <- ncol(g0covs)-1
      recovs <- reForm$recovs
      nbRecovs <- ncol(recovs)-1
      if(missing(covs)){
        covs <- model$data
      } else if(is.matrix(covs)){
        stop("covs must be provided as a data frame for recharge models")
      }

      if(inherits(model,"hierarchical")) class(covs) <- append("hierarchical",class(covs))
      tmpSplineInputs<-getSplineFormula(newformula,model$data,covs)
      
      testCovs <- tmpSplineInputs$covs
      if(inherits(model,"hierarchical")){
        testCovs$level <- model$data$level[1]
      }
      # check that all covariates are provided
      ck1 <- tryCatch(stats::model.matrix(recharge$theta,testCovs),error=function(e) e)
      ck2 <- tryCatch(stats::model.matrix(tmpSplineInputs$formula,testCovs),error=function(e) e)
      if(inherits(ck1,"error")) stop("covs not specified correctly -- ",ck1)
      if(inherits(ck2,"error")) stop("covs not specified correctly -- ",ck2)
    } else {
      if(missing(covs)){
        covs <- model$rawCovs
        if(length(covs)){
          tmpSplineInputs<-getSplineFormula(newformula,model$data,covs)
        } else {
          tmpSplineInputs<-getSplineFormula(newformula,model$data,model$data[1,])
        }
      }
    }
    
    if(!is.matrix(covs)){
      if(inherits(model,"hierarchical") & !("level" %in% names(covs))){
        # expand covs for each level of hierarchy
        tmpSplineInputs$covs <- data.frame(tmpSplineInputs$covs[rep(1:nrow(tmpSplineInputs$covs),nlevels(model$data$level)),,drop=FALSE],level=rep(levels(model$data$level),each=nrow(tmpSplineInputs$covs)))
        class(tmpSplineInputs$covs) <- append("hierarchical",class(tmpSplineInputs$covs))
      }
      covMat <- stats::model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
      if(!is.null(covIndex)) {
        if(!is.numeric(covIndex) || any(covIndex<1 | covIndex>nrow(covMat))) stop("covIndex can only include integers between 1 and ",nrow(covMat))
      } else covIndex <- 1:nrow(covMat)
    }
    
    mixtures <- model$conditions$mixtures
  
    probs <- list()
      
    for(mix in 1:mixtures){
      # all transition matrices
      tmpbeta <- beta[(mix-1)*ncol(covMat)+1:ncol(covMat),,drop=FALSE]
      if(is.null(recharge)){
        if(!is.null(covIndex)) {
          if(!is.numeric(covIndex) || any(covIndex<1 | covIndex>nrow(covMat))) stop("covIndex can only include integers between 1 and ",nrow(covMat))
        } else covIndex <- 1:nrow(covMat)
        allMat <- trMatrix_rcpp(nbStates=nbStates, beta=tmpbeta, covs=covMat[covIndex,,drop=FALSE], betaRef=model$conditions$betaRef)
      } else {
        if(!is.null(covIndex)) {
          if(!is.numeric(covIndex) || any(covIndex<1 | covIndex>nrow(tmpSplineInputs$covs))) stop("covIndex can only include integers between 1 and ",nrow(tmpSplineInputs$covs))
        } else covIndex <- 1:nrow(tmpSplineInputs$covs)
        gamInd<-(length(model$mod$estimate)-(nrow(tmpbeta))*nbStates*(nbStates-1)*mixtures+1):(length(model$mod$estimate))-(ncol(model$covsPi)*(mixtures-1))-ifelse(nbRecovs,(nbRecovs+1)+(nbG0covs+1),0)-ncol(model$covsDelta)*(nbStates-1)*(!model$conditions$stationary)*mixtures
        allMat <- array(unlist(lapply(split(tmpSplineInputs$covs[covIndex,,drop=FALSE],covIndex),function(x) tryCatch(get_gamma_recharge(model$mod$estimate[c(gamInd[unique(c(model$conditions$betaCons))],length(model$mod$estimate)-nbRecovs:0)],covs=x,formula=tmpSplineInputs$formula,hierRecharge=hierRecharge,nbStates=nbStates,betaRef=model$conditions$betaRef,betaCons=model$conditions$betaCons,workBounds=rbind(model$conditions$workBounds$beta,model$conditions$workBounds$theta),mixture=mix),error=function(e) NA))),dim=c(nbStates,nbStates,nrow(tmpSplineInputs$covs)))
      }
      
      tryCatch({
        
          # for each transition matrix, derive corresponding stationary distribution
          if(!inherits(model,"hierarchical")){

            probs[[mix]] <- getProbs(allMat,model$stateNames)
          
          } else {
            
            probs[[mix]] <- list()
            for(j in 1:(model$conditions$hierStates$height-1)){
                  
              if(j==1){
                ref <- model$conditions$hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)
                probs[[mix]][["level1"]] <- getProbs(allMat[ref,ref,which(covMat[covIndex,colnames(covMat) %in% paste0("I((level == \"",j,"\") * 1)")]==1),drop=FALSE],names(ref))
              } else {
                
                t <- data.tree::Traverse(model$conditions$hierStates,filterFun=function(x) x$level==j)
                names(t) <- model$conditions$hierStates$Get("name",filterFun=function(x) x$level==j)
                
                if(length(names(t))) probs[[mix]][[paste0("level",j)]] <- list()
                
                for(k in names(t)){
                  ref <- t[[k]]$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)#t[[k]]$Get("state",filterFun = data.tree::isLeaf)
                  if(!is.null(ref)){
                    probs[[mix]][[paste0("level",j)]][[k]] <- getProbs(allMat[ref,ref,which(covMat[covIndex,colnames(covMat) %in% paste0("I((level == \"",j,"\") * 1)")]==1),drop=FALSE],names(ref))
                  }
                }  
              }
            }
          }
      },
      error = function(e) {
          stop(paste("The stationary probabilities cannot be calculated",
                     "for these covariate values (singular system)."))
      })
    }
    return(probs)
}

#' @method stationary miSum
#' @export
stationary.miSum <- function(model, covs, covIndex=NULL)
{
  model <- formatmiSum(model)

  stationary(momentuHMM(model), covs, covIndex)
}

#' @method stationary miHMM
#' @export
stationary.miHMM <- function(model, covs, covIndex=NULL)
{
  stationary(model$miSum, covs, covIndex)
}

getProbs <- function(allMat,stateNames){
  
  nbStates <- length(stateNames)
  
  probs <- apply(allMat, 3,
                 function(gamma)
                     solve(t(diag(nbStates)-gamma+1),rep(1,nbStates)))
  probs <- t(probs)
  
  colnames(probs) <- stateNames
  
  probs
}

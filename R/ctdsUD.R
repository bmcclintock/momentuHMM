#' Calculate and plot the (state-dependent) stationary (or utilization) distribution from a \code{ctds} fitted model object
#' 
#' This is essentially a wrapper function for \code{\link[ctmcmove]{get.rate.matrix}} and \code{\link[ctmcmove]{get.UD}} that has been extended for \code{\link{momentuHMM}}. 
#' 
#' @param ctds A \code{ctds} object returned by \code{\link{fitCTHMM}}
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatially referenced covariates. Raster \code{\link[raster]{stack}}s or \code{\link[raster]{brick}}s are not currently allowed.
#' @param spatialCovs.grad List of \code{\link[raster]{raster}} objects for spatially referenced covariates, where a directional gradient is to be calculated internally using \code{\link[ctmcmove]{rast.grad}}. Raster \code{\link[raster]{stack}}s or \code{\link[raster]{brick}}s are not currently allowed. 
#' @param zero.idx Integer vector of the indices of raster cells that are not passable and should be excluded. These are cells where movement should be impossible. Default is zero.idx=integer().
#' @param method	Either "lu" (default) or "limit". See \code{\link[ctmcmove]{get.UD}} details for a description of the two methods.
#' @param maxiter	Total number of iterations for limit method if tolerance not reached first. Defaults to 100. Ignored for method='lu'.
#' @param start	A value for the starting distribution for the 'limit' method. Defaults to 1/num. cells. Ignored for method='lu'.
#' @param tol	Value used to assess convergence for limit method. If max(abs(pi1-pi0))<tol, limit method has converged. Defaults to sqrt(.Machine$double.eps).
#' @return A list of \code{\link[raster]{raster}} objects.
#' @export
ctdsUD <- function (ctds, spatialCovs, spatialCovs.grad, zero.idx = integer(), method="lu",maxiter, start, tol) 
{
  if(inherits(ctds,"miHMM")) ctds <- ctds$miSum
  if(inherits(ctds,"miSum")){
    ctds <- formatmiSum(ctds)
    ctds$CIbeta <- ctds$Par$beta
  }
  
  if(!inherits(ctds,"ctds")) stop("ctds must be a 'ctds' object")
  
  if(missing(spatialCovs)) spatialCovs <- NULL
  if(missing(spatialCovs.grad)) spatialCovs.grad <- NULL
  
  if(is.null(spatialCovs) & is.null(spatialCovs.grad)) stop("No spatial covariates provided!")
  checkRast(ctds$data[,which(!names(ctds$data) %in% c(names(spatialCovs),names(spatialCovs.grad)))],spatialCovs,spatialCovs.grad)
  
  directions <- attr(ctds$data,"directions")
  normalize.gradients <- attr(ctds$data,"normalize.gradients")
  grad.point.decreasing <- attr(ctds$data,"grad.point.decreasing")
  
  if(length(spatialCovs)){
    examplerast <- spatialCovs[[1]]
  } else examplerast <- spatialCovs.grad[[1]]
  
  nn = raster::ncell(examplerast)
  locs = 1:nn
  adj = raster::adjacent(examplerast, locs, pairs = TRUE, sorted = TRUE, 
                 id = TRUE, directions = directions)
  idx.mot = adj[, 2:3]
  
  start.cells = idx.mot[, 1]
  adj.cells = idx.mot[, 2]
  xy.cell = raster::xyFromCell(examplerast, start.cells)
  xy.adj = raster::xyFromCell(examplerast, adj.cells)
  v.adj = (xy.adj - xy.cell)/sqrt(apply((xy.cell - xy.adj)^2, 1, sum))
  
  X <- NULL
  if (length(spatialCovs)) {
    if(any(unlist(lapply(spatialCovs,raster::nlayers))>1)){
      stop("spatialCovs elements cannot contain more than one layer; if this is a spatio-temporal covariate, only one of these layer can be supplied at time")
    } else {
      X.static <- get.static(directions=directions,start.cells=start.cells,spatialCovs=spatialCovs)
    }
    X <- data.frame(X.static)
  } else p.static <- 0 
  
  p.crw = 0
  
  if (length(spatialCovs.grad)) {
    if(any(unlist(lapply(spatialCovs.grad,raster::nlayers))>1)){
      stop("spatialCovs.grad elements cannot contain more than one layer; if this is a spatio-temporal covariate, only one of these layer can be supplied at time")
    } else {
      X.grad <- get.grad(directions=directions,start.cells=start.cells,v.adj=v.adj,spatialCovs.grad=spatialCovs.grad,normalize.gradients=normalize.gradients,grad.point.decreasing =  grad.point.decreasing)
    }
    if(is.null(X)) X <- data.frame(X.grad)
    else X <- cbind(X,X.grad)
  } else  p.grad = 0
  
  X$tau = 1
  X$crw = 0
  
  
  stateNames <- ctds$stateNames
  nbStates <- length(stateNames)
  R <- udrast <- vector('list',nbStates)
  names(R) <- names(udrast) <- stateNames
  form <- stateFormulas(ctds$conditions$DM[[attr(ctds$data,"ctdsData")]]$lambda,nbStates)
  #parInd <- apply(matrix(unlist(lapply(ctds$conditions$fullDM$z,function(x) all(x==0))),ncol=ncol(ctds$conditions$fullDM$z))[1:nbStates,,drop=FALSE],1,function(y) which(y==0))
  ##if(nbStates==1) parInd <- list(parInd)
  for(i in 1:nbStates){
    R[[stateNames[[i]]]] = Matrix::Matrix(0, nrow = nn, ncol = nn, sparse = TRUE)
    DM <- stats::model.matrix(form[[i]],X)
    parInd <- which(!unlist(lapply(ctds$conditions$fullDM$z[i,],function(x) all(x==0))))
    if("crw" %in% all.vars(form[[i]])) {
      warning("It is not possible to include 'crw' when calculating UD for state ",i,"; this term has been set to zero")
      DM <- DM[,-which(colnames(DM)=="crw"),drop=FALSE]
      if(ncol(DM)==1) next;
    }
    R[[stateNames[i]]][idx.mot] <- exp(DM %*% ctds$CIbeta[[attr(ctds$data,"ctdsData")]]$est[parInd])
    if (length(zero.idx) > 0) {
      R[[stateNames[i]]][zero.idx, zero.idx] = 0
    }
    R[[stateNames[i]]] <- ctmcmove::get.UD(R[[stateNames[i]]], method=method,maxiter=maxiter, start=start, tol=tol)
    udrast[[stateNames[i]]] <- examplerast
    raster::values(udrast[[stateNames[i]]]) <- as.numeric(R[[stateNames[i]]])
    names(udrast[[stateNames[i]]]) <- "UD"
    raster::plot(udrast[[stateNames[i]]],main=stateNames[i])
  }
  return(udrast)
}

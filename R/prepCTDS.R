
#' Preprocessing of continuous-time discrete-space (CTDS) movement HMMs using ctmcmove
#' 
#' This wrapper function for \code{\link[ctmcmove]{path2ctmc}} and \code{\link[ctmcmove]{ctmc2glm}} converts a \code{data.frame} of coordinates, other data streams, and non-spatial covariates to a \code{\link{momentuHMMData}} object that can be passed directly to \code{\link{fitCTMM}} (or as a list to \code{\link{MIfitCTHMM}}).
#' 
#' @param data A \code{data.frame} that must include entries for the x coordinate (\code{x}), y coordinate (\code{y}), and time stamp (\code{time}). An \code{ID} entry must also be included if \code{data} includes multiple individuals.
#' @param rast A raster object or raster stack object that will define the discrete-space grid cells for the CTMC movement path.
#' @param directions Integer. Either 4 (indicating a "Rook's neighborhood" of 4 neighboring grid cells) or 8 (indicating a "King's neighborhood" of 8 neighboring grid cells).
#' @param zero.idx Integer vector of the indices of raster cells that are not passable and should be excluded. These are cells where movement should be impossible. Default is zero.idx=integer().
#' @param print.iter Logical. If true, then the progress stepping through each observed location in "xy" and "t" will be output in the terminal.
#' @param method Specifies interpolation method. Either "ShortestPath", which uses the shortest graphical path on the raster graph, or "LinearInterp", which linearly interpolates between observed locations. "ShortestPath" is slower, slightly more accurate, and allows for impassible barriers specified through "zero.idx". "LinearInterp" is faster but does not allow for impassible barriers.
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'Date'). In the \code{\link{momentuHMMData}} object returned by \code{prepCTDS}, covariates for the current position (e.g.\ for use in \code{formula} or \code{DM}) are named with a \code{.cur} suffix (e.g. \code{cov1.cur}).
#' @param spatialCovs.grad List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates, where a directional gradient is to be calculated internally using \code{\link[ctmcmove]{rast.grad}}. Gradient-based covariates specified by \code{spatialCovs.grad} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs.grad} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'Date').
#' @param crw	Logical. If TRUE (default), an autocovariate is created for each cell visited in the CTMC movement path. The autocovariate is a unit-length directional vector pointing from the center of the most recent grid cell to the center of the current grid cell.
#' @param normalize.gradients	Logical. Default is FALSE. If TRUE, then all gradient covariates for \code{spatialCovs.grad} are normalized by dividing by the length of the gradient vector at each point.
#' @param grad.point.decreasing	Logical. If TRUE, then the gradient covariates are positive in the direction of decreasing values of the covariate. If FALSE, then the gradient covariates are positive in the direction of increasing values of the covariate (like a true gradient).
#' @param covNames Character vector indicating the names of any covariates in \code{data}. Any variables in \code{data} (other than \code{ID}, \code{x}, \code{y}, and \code{time}) that are not identified in covNames are assumed to be additional data streams (i.e., missing values will not be accounted for).
#' @return A \code{\link{momentuHMMData}} object
#' @export
prepCTDS <- function(data, rast, directions=4, zero.idx=integer(), print.iter=FALSE, method="ShortestPath",
                     spatialCovs=NULL, spatialCovs.grad=NULL, crw = TRUE, normalize.gradients = FALSE, grad.point.decreasing = TRUE,
                     covNames=NULL) {
  
  if (!requireNamespace("ctmcmove", quietly = TRUE)) {
    stop("Package \"ctmcmove\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if(is.null(data$x) | is.null(data$y) | is.null(data$time))
    stop("data must contain 'x', 'y', and 'time' fields")
  
  if(any(names(data) %in% c("x.current","y.current","z","itTPM"))) stop("'x.current', 'y.current', 'z', and 'itTPM' are reserved and cannot be fields in data")
  if(!is.null(covNames)){
    if(any(covNames %in% c(c("ID","time","x","y")))) stop("covNames cannot include 'ID', 'time', 'x', or 'y'")
  }
  
  data <- prepData(data,coordNames=NULL,covNames=covNames)
  dataStreams <- names(data)[which(!names(data) %in% c("ID","time","x","y",covNames))]
  
  # check rasters using prepData
  check1 <- prepData(data[1:3,],spatialCovs=spatialCovs)
  check2 <- prepData(data[1:3,],spatialCovs=spatialCovs.grad)
  if(any(names(spatialCovs) %in% names(spatialCovs.grad))) stop("'spatialCovs' and 'spatialCovs.grad' names must be unique")
  
  ctdsglm <- list()
  for(i in unique(data$ID)){
    ctds <- path2ctds(xy=as.matrix(data[which(data$ID==i & !is.na(data$x) & !is.na(data$y)),c("x","y")]),t=data$time[which(data$ID==i & !is.na(data$x) & !is.na(data$y))],rast=rast,directions=directions,zero.idx=zero.idx,print.iter=print.iter,method=method)
    ctdsglm[[i]] <- ctds2glm(data[which(data$ID==i & !is.na(data$x) & !is.na(data$y)),], ctds,rast = rast, directions=directions, spatialCovs = spatialCovs, spatialCovs.grad=spatialCovs.grad, crw=crw, normalize.gradients = normalize.gradients, grad.point.decreasing = grad.point.decreasing, include.cell.locations = TRUE, zero.idx=zero.idx, covNames = covNames)
    ## add back initial time so dt can be correctly calculated in fitCTHMM (and covariates can be drawn based on the initial position)
    #tmp <- ctdsglm[[i]][1:directions,]
    #tmp$t <- data[which(data$ID==i),"time"][1]
    #tmp$z <- NA
    #tmp[,c("x.current","y.current")] <- rep(xyFromCell(int,cellFromXY(int,data[which(data$ID==i)[1],c("x","y")])),each=directions)
    #ctdsglm[[i]] <- rbind(tmp,ctdsglm[[i]])
    if(length(dataStreams)){
      multiCellMove <- which(ctdsglm[[i]]$cellCross>0)
      if(length(multiCellMove)){
        warning("There are moves across multiple cells within a time step: \n",
        "   -- any spatial covariates for ",paste0(dataStreams,collapse=", ")," pertain to the initial cell for these time step(s) \n",
        "   -- during model fitting, the state(s) for ",paste0(dataStreams,collapse=", ")," are assumed to be the initial state(s) at the start of these time step(s)")
        for(j in unique(ctdsglm[[i]]$cellCross[multiCellMove])){
          crInd <- which(ctdsglm[[i]]$cellCross==j)
          tmp <- ctdsglm[[i]][crInd,][1:directions,]
          tmp$tau <- sum(ctdsglm[[i]][which(ctdsglm[[i]]$cellCross==j),"tau"])/directions
          tmp$z <- NA
          ctdsglm[[i]][which(ctdsglm[[i]]$cellCross==j),dataStreams] <- NA
          ctdsglm[[i]] <- rbind(ctdsglm[[i]][1:(crInd[1]-1),],tmp,ctdsglm[[i]][-(1:(crInd[1]-1)),])
        }
      }
      ctdsglm[[i]][(1:nrow(ctdsglm[[i]]))[-seq(1,nrow(ctdsglm[[i]]),directions)],dataStreams] <- NA
    }
    names(ctdsglm[[i]])[which(names(ctdsglm[[i]])=="t")] <- "time"
  }
  ctdsglm <- do.call(rbind,ctdsglm)
  rownames(ctdsglm) <- NULL
  
  ctdsglm <- ctdsglm[,which(!colnames(ctdsglm) %in% c("x.adj","y.adj"))]
  #ctdsOut <- ctdsglm[seq(1,nrow(ctdsglm),directions),]
  #names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x.current","y.current",dataStreams,covNames))] <- paste0(names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x.current","y.current",dataStreams,covNames))],".1")
  #for(j in 2:directions){
  #  tmp <- ctdsglm[seq(j,nrow(ctdsglm),directions),which(!colnames(ctdsglm) %in% c("ID","step","angle","time","x.current","y.current",dataStreams,covNames))]
  #  names(tmp) <- paste0(names(tmp),".",j)
  #  ctdsOut <- cbind(ctdsOut,tmp)
  #}
  
  # add non-gradient spatial covariates for current position (e.g. for inclusion in TPM)
  #ctdsOut <- prepData(ctdsOut,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs,altCoordNames = "current")
  names(spatialCovs) <- paste0(names(spatialCovs),".cur")
  ctdsOut <- prepData(ctdsglm,coordNames=c("x.current","y.current"),spatialCovs=spatialCovs)
  ctdsOut <- ctdsOut[,c("ID","time","x","y","z",dataStreams,covNames,names(ctdsOut)[which(!names(ctdsOut) %in% c("ID","time","x","y","step","angle","z",dataStreams,covNames))])]
  
  class(ctdsOut) <- unique(append(c("momentuHMMData","ctds"),class(ctdsOut)))
  attr(ctdsOut,"directions") <- directions
  attr(ctdsOut,"coords") <- c("x","y")
  attr(ctdsOut,"prodPois") <- "z"
  return(ctdsOut)
}

path2ctds<-function (xy, t, rast, directions = 4, zero.idx = integer(), 
          print.iter = FALSE, method = "ShortestPath") 
{
  if (class(rast) == "RasterStack") {
    rast = rast[[1]]
  }
  raster::values(rast) <- 1
  raster::values(rast)[zero.idx] <- 0
  trans = gdistance::transition(rast, prod, directions = directions)
  ncell = raster::ncell(rast)
  if (method == "LinearInterp") {
    A = Matrix::Matrix(0, nrow = ncell, ncol = ncell, sparse = TRUE)
    adj = raster::adjacent(rast, 1:ncell)
    A[adj] <- 1
  }
  path = cbind(xy, t)
  tidx = sort(t, index.return = T)$ix
  path = path[tidx, ]
  T = nrow(path)
  ec.all = raster::cellFromXY(rast, xy)
  ec = ec.all[1]
  current.cell = ec
  rt = integer()
  current.rt = 0
  cellCross <- integer()
  current.cr <- 0
  if (print.iter) {
    cat("Total locations =", T, "\n")
  }
  for (i in 2:(T)) {
    if (print.iter) {
      cat(i, " ")
    }
    if (ec.all[i] == current.cell) {
      current.rt = (t[i] - t[i - 1])
      rt <- c(rt,(t[i] - t[i - 1]))
      ec = c(ec, ec.all[i])
      cellCross <- c(cellCross,0)
    }
    else {
      if (method == "ShortestPath") {
        sl = gdistance::shortestPath(trans, as.numeric(xy[i - 1, 
        ]), as.numeric(xy[i, ]), "SpatialLines")
        slc = raster::coordinates(sl)[[1]][[1]]
        sl.cells = raster::cellFromXY(rast, slc)
        t.in.each.cell = (t[i] - t[i - 1])/length(sl.cells)
        rt = c(rt, rep(t.in.each.cell, length(sl.cells)))
        current.rt = t.in.each.cell
        ec = c(ec, sl.cells)
        current.cell = ec[length(ec)]
        current.cr <- current.cr + 1
        cellCross <- c(cellCross,rep(current.cr,length(sl.cells)))
      }
      if (method == "LinearInterp") {
        if (A[current.cell, ec.all[i]] == 1) {
          rt = c(rt, (t[i] - t[i - 1]))
          current.rt = 0
          ec = c(ec, ec.all[i])
          current.cell = ec.all[i]
          cellCross <- c(cellCross,0)
        }
        else {
          xyt.1 = path[i - 1, ]
          xyt.2 = path[i, ]
          d = sqrt(sum((xyt.1[-3] - xyt.2[-3])^2))
          rast.res = raster::res(rast)[1]
          xapprox = stats::approx(c(xyt.1[3], xyt.2[3]), c(xyt.1[1], 
                                                    xyt.2[1]), n = max(100, round(d/rast.res * 
                                                                                    100)))
          yapprox = stats::approx(c(xyt.1[3], xyt.2[3]), c(xyt.1[2], 
                                                    xyt.2[2]), n = max(100, round(d/rast.res * 
                                                                                    100)))
          tapprox = xapprox$x
          xapprox = xapprox$y
          yapprox = yapprox$y
          xycell.approx = raster::cellFromXY(rast, cbind(xapprox, 
                                                 yapprox))
          rle.out = rle(xycell.approx)
          ec.approx = rle.out$values
          transition.idx.approx = c(1,cumsum(rle.out$lengths))
          transition.times.approx = tapprox[transition.idx.approx]
          rt.approx = diff(transition.times.approx)#diff(c(t[i-1],transition.times.approx))#
          rt = c(rt, rt.approx)
          current.rt = rt.approx[length(rt.approx)]
          ec = c(ec, ec.approx[-1],ec.approx[length(ec.approx)])
          current.cell = ec[length(ec)]
          current.cr <- current.cr + 1
          cellCross <- c(cellCross,rep(current.cr,length(rt.approx)))
        }
      }
    }
  }
  list(ec = ec, rt = c(diff(c(t[1],t[1] + cumsum(rt))),0), trans.times = c(t[1],t[1] + cumsum(rt)), cellCross=c(cellCross,0))
}

ctds2glm <- function (data, ctmc, rast, spatialCovs=NULL, spatialCovs.grad=NULL, crw = TRUE, normalize.gradients = FALSE, 
          grad.point.decreasing = TRUE, include.cell.locations = TRUE, 
          directions = 4, zero.idx = integer(), covNames) 
{
  if(length(spatialCovs)){
    examplerast <- spatialCovs[[1]]
  } else examplerast <- rast
  
  locs = ctmc$ec
  wait.times = ctmc$rt
  notzero.idx = 1:raster::ncell(examplerast)
  if (length(zero.idx) > 0) {
    notzero.idx = notzero.idx[-zero.idx]
  }
  adj = raster::adjacent(examplerast, locs, pairs = TRUE, sorted = TRUE, 
                         id = TRUE, directions = directions, target = notzero.idx)
  adj.cells = adj[, 3]
  rr = rle(adj[, 1])
  time.idx = rep(rr$values, times = rr$lengths)
  start.cells = adj[, 2]
  moves <- which(locs[-length(locs)]!=locs[-1])
  z = rep(0, length(start.cells))
  idx.move = rep(0, length(z))
  diag.move = rep(0, length(locs))
  for (i in moves) {
    idx.t = which(time.idx == i)
    idx.m = which(adj.cells[idx.t] == locs[i + 1])
    z[idx.t[idx.m]] <- 1
    if (length(idx.m) == 0) {
      diag.move[i] = 1
      #z[idx.t] <- NA
      error("diagonal move detected and cannot be properly accounted for",ifelse(directions==4,"; try expanding 'directions' to 8",""))
    }
  }
  tau = rep(wait.times, times = rr$lengths)
  t = rep(ctmc$trans.times, times = rr$lengths)
  cr <- rep(ctmc$cellCross, times = rr$length)
  
  xy.cell = raster::xyFromCell(examplerast, start.cells)
  xy.adj = raster::xyFromCell(examplerast, adj.cells)
  v.adj = (xy.adj - xy.cell)/sqrt(apply((xy.cell - xy.adj)^2, 
                                        1, sum))
  
  moveData <- merge(data.frame(ID=data$ID[1],time=ctmc$trans.times),data,by=c("ID","time"),all=TRUE)
  moveData <- suppressWarnings(prepData(moveData,coordNames=NULL,covNames=covNames))
  moveData <- moveData[which(moveData$time %in% ctmc$trans.times),]
  dataStreams <- names(moveData)[which(!names(moveData) %in% c("ID","time","x","y",covNames))]
  
  ### extract spatial covariates
  # gradient-based covariates
  p.grad = length(spatialCovs.grad)

  if(p.grad) {
    X.grad = do.call(cbind,lapply(spatialCovs.grad,
                function(x){
                  if(inherits(x,"RasterLayer")){
                    x = ctmcmove::rast.grad(x)
                    if (normalize.gradients) {
                      lengths = sqrt(x$grad.x^2 + x$grad.y^2)
                      x$grad.x <- x$grad.x/lengths
                      x$grad.y <- x$grad.y/lengths
                    }
                    return(v.adj[, 1] * x$grad.x[start.cells] + v.adj[, 2] * x$grad.y[start.cells])
                  } else {
                    zname <- names(attributes(x)$z)
                    zvalues <- raster::getZ(x)
                    if(inherits(x,"RasterBrick")) x <- raster::setZ(raster::stack(x),zvalues,zname)
                    grad <- ctmcmove::rast.grad(x[[1]])
                    gradx <- raster::stack(grad$rast.grad.x)
                    grady <- raster::stack(grad$rast.grad.y)
                    for(i in 2:nlayers(x)){
                      grad <- ctmcmove::rast.grad(x[[i]])
                      gradx <- raster::stack(gradx,grad$rast.grad.x)
                      grady <- raster::stack(grady,grad$rast.grad.y)
                    }
                    names(gradx) <- names(grady) <- names(x)
                    gradx <- raster::setZ(gradx,zvalues,zname)
                    grady <- raster::setZ(grady,zvalues,zname)
                    fullx <- gradx[start.cells]
                    fully <- grady[start.cells]
                    tmpspCovs.x <- tmpspCovs.y <- numeric(length(start.cells))
                    if((!zname %in% covNames) & zname!="time") stop("z-value name '",zname,"' must be included in covNames")
                    for(ii in 1:length(zvalues)){
                      tmpspCovs.x[which(rep(moveData[[zname]],each=directions)==zvalues[ii])] <- fullx[which(rep(moveData[[zname]],each=directions)==zvalues[ii]),ii]
                      tmpspCovs.y[which(rep(moveData[[zname]],each=directions)==zvalues[ii])] <- fully[which(rep(moveData[[zname]],each=directions)==zvalues[ii]),ii]
                    }
                    return(v.adj[, 1] * tmpspCovs.x + v.adj[, 2] * tmpspCovs.y)
                  }
                }))
    
    if (grad.point.decreasing == TRUE) {
      X.grad = -X.grad
    }
  }
  
  # non-gradient based covariates
  X.static = do.call(cbind,lapply(spatialCovs,
                function(x){
                  fullspCovs <- x[start.cells]
                  if(inherits(x,"RasterLayer")){
                    return(fullspCovs)
                  } else {
                    spCov <- numeric(length(start.cells))
                    zname <- names(attributes(x)$z)
                    if((!zname %in% covNames) & zname!="time") stop("z-value name '",zname,"' must be included in covNames")
                    zvalues <- raster::getZ(x)
                    for(ii in 1:length(zvalues)){
                      spCov[which(rep(moveData[[zname]],each=directions)==zvalues[ii])] <- fullspCovs[which(rep(moveData[[zname]],each=directions)==zvalues[ii]),ii]
                    }
                    return(spCov)
                  }
                }))
  
  zInd <- which(rowSums(matrix(z,ncol=directions,byrow=TRUE))==1)
  if(zInd[length(zInd)]*directions<length(z)) zInd <- c(zInd,zInd[length(zInd)]+1)
  moveInd = c(t(matrix(1:length(z),ncol=directions,byrow=TRUE)[zInd,]))
  X.crw <- numeric(length(z))
  idx.move = which(z == 1)
  idx.move = c(idx.move, min(length(z),tail(idx.move,1)+directions))
  v.moves = v.adj[rep(idx.move[1:(length(rr$lengths[zInd]) - 1)], 
                      times = rr$lengths[zInd][-1]), ]
  v.moves = rbind(matrix(0, ncol = 2, nrow = directions), 
                  v.moves)
  X.crw[moveInd] = apply(v.moves * v.adj[moveInd,], 1, sum)
  zInd <- c(0,zInd)
  for(i in 1:(length(zInd)-1)){
    if(zInd[i+1]-zInd[i] > 1) {
      X.crw[(zInd[i]*directions+1):((zInd[i+1]-1)*directions)] <- rep(matrix(X.crw,ncol=directions,byrow=TRUE)[zInd[i+1],],(zInd[i+1]-1)-zInd[i])
    }
  }
  X.crw[((zInd[length(zInd)])*directions+1):length(X.crw)] <- matrix(X.crw,ncol=directions,byrow=TRUE)[zInd[i+1],]
  
  if (crw == FALSE & p.grad > 0) {
    X = cbind(X.static, X.grad)
  }
  if (crw == TRUE & p.grad > 0) {
    X = cbind(X.static, X.grad, X.crw)
    colnames(X)[ncol(X)] = "crw"
  }
  if (crw == FALSE & p.grad == 0) {
    X = cbind(X.static)
  }
  if (crw == TRUE & p.grad == 0) {
    X = cbind(X.static, X.crw)
    colnames(X)[ncol(X)] = "crw"
  }
  if (include.cell.locations) {
    xys = cbind(xy.cell, xy.adj)
    colnames(xys) = c("x.current", "y.current", "x.adj", 
                      "y.adj")
    X = cbind(X, xys)
  }
  if(length(dataStreams)) X = cbind(X,moveData[,dataStreams,drop=FALSE][rep(seq_len(nrow(moveData)),each=directions),,drop=FALSE])
  if(length(covNames)) X = cbind(X,moveData[,covNames,drop=FALSE][rep(seq_len(nrow(moveData)),each=directions),,drop=FALSE])
  out = data.frame(z = z, X, tau = tau, t = t, cellCross = cr)
  T = nrow(out)
  out = out[-((T - (directions-1)):T), ]
  out
}

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

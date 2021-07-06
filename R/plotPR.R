
#' Plot pseudo-residuals
#'
#' Plots time series, qq-plots (against the standard normal distribution) using \code{\link[car]{qqPlot}}, and sample
#' ACF functions of the pseudo-residuals for each data stream
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{momentuHierHMM}}, \code{\link{miHMM}}, \code{\link{HMMfits}}, or \code{\link{miSum}} object.
#' @param lag.max maximum lag at which to calculate the acf.  See \code{\link[stats]{acf}}.
#' @param ncores number of cores to use for parallel processing
#'
#' @details \itemize{
#' \item If some turning angles in the data are equal to pi, the corresponding pseudo-residuals
#' will not be included. Indeed, given that the turning angles are defined on (-pi,pi], an angle of pi
#' results in a pseudo-residual of +Inf (check Section 6.2 of reference for more information on the
#' computation of pseudo-residuals).
#' \item If some data streams are zero-inflated and/or one-inflated, the corresponding pseudo-
#' residuals are shown as segments, because pseudo-residuals for discrete data are defined as
#' segments (see Zucchini and MacDonald, 2009, Section 6.2).
#' \item For multiple imputation analyses, if \code{m} is a \code{\link{miHMM}} object or a list of \code{\link{momentuHMM}} objects, then
#' the pseudo-residuals are individually calculated and plotted for each model fit. Note that pseudo-residuals for \code{\link{miSum}} objects (as returned by \code{\link{MIpool}}) are based on pooled parameter 
#' estimates and the means of the data values across all imputations (and therefore may not be particularly meaningful).
#' }
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plotPR(m)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export
#' @importFrom stats acf na.pass qqnorm
# #' @importFrom car qqPlot

plotPR <- function(m, lag.max = NULL, ncores = 1)
{
  
  m <- delta_bc(m)
  
  nSims <- 1
  if(!is.momentuHMM(m) & !is.miSum(m)){
    if(!is.miHMM(m) & !is.HMMfits(m)) stop("'m' must be a momentuHMM, HMMfits, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
    else {
      if(is.miHMM(m)) m <- m$HMMfits
      goodIndex <- which(unlist(lapply(m,is.momentuHMM)))
      if(length(goodIndex) < length(m))  warning("The following imputations are not momentuHMM objects and will be ignored: ",paste((1:length(m))[-goodIndex],collapse=", "))
      m <- m[goodIndex]
      distnames <- names(m[[1]]$conditions$dist)
      nSims <- length(m)
    }
  } else distnames <- names(m$conditions$dist)

  cat("Computing pseudo-residuals",ifelse(nSims>1,paste0(" for ",nSims," model fits... \n"),"... "),sep="")
  pr <- pseudoRes(m, ncores = ncores)
  
  if(nSims==1){
    pr <- list(pr)
    m <- list(m)
    col <- "red"
  } else {
    hues <- seq(15, 375, length = nSims + 1)
    col <- c("red",hcl(h = hues, l = 65, c = 100)[1:nSims])
  }
  par(mfrow=c(length(distnames),3))

  # reduce top margin
  par(mar=c(5,4,4,2)-c(0,0,3,0)) # bottom, left, top, right
  
  if(any(unlist(lapply(m,function(x) inherits(x,"hierarchical"))))){
    installDataTree()
    for(i in distnames){
      for(j in 1:nSims){
        iLevel <- gsub("level","",m[[j]]$conditions$hierDist$Get("parent",filterFun=data.tree::isLeaf)[[i]]$name)
        pr[[j]][[paste0(i,"Res")]] <- pr[[j]][[paste0(i,"Res")]][which(m[[j]]$data$level==iLevel)]
      }
    }
  }
  
  for(i in distnames){
      # time series
      if(any(is.infinite(pr[[1]][[paste0(i,"Res")]]))){
        warning("some pseudo-residuals for ",i," are infinite; these will be set to NA for plotting")
        pr[[1]][[paste0(i,"Res")]][which(is.infinite(pr[[1]][[paste0(i,"Res")]]))] <- NA
      }
      if(all(is.na(pr[[1]][[paste0(i,"Res")]]))) next;
      ylim <- range(unlist(lapply(pr,function(x) x[[paste0(i,"Res")]][which(is.finite(x[[paste0(i,"Res")]]))])),na.rm=TRUE)
      plot(pr[[1]][[paste0(i,"Res")]],type="l",xlab="Observation index",ylab=paste0(i," pseudo-residuals"),main="",ylim=ylim)
      if(nSims>1){
        for(j in 2:nSims){
          if(any(is.infinite(pr[[j]][[paste0(i,"Res")]]))){
            warning("some pseudo-residuals for ",i," are infinite; these will be set to NA for plotting")
            pr[[j]][[paste0(i,"Res")]][which(is.infinite(pr[[j]][[paste0(i,"Res")]]))] <- NA
          }
          lines(pr[[j]][[paste0(i,"Res")]],col=col[j])
        }
      }
      
      # qq-plot
      if(m[[1]]$conditions$dist[[i]] %in% mvndists){
        tpr <- as.matrix(pr[[1]][[paste0(i,"Res")]])
        if(m[[1]]$conditions$dist[[i]]=="mvnorm2" || m[[1]]$conditions$dist[[i]]=="mvnorm3") ndim <- as.numeric(gsub("mvnorm","",m[[1]]$conditions$dist[[i]] ))
        else ndim <- as.numeric(gsub("rw_mvnorm","",m[[1]]$conditions$dist[[i]] ))
        if (!requireNamespace("car", quietly = TRUE)) {
          stop("Package \"car\" needed for QQ plot of mutivariate normal distributions. Please install it.",
               call. = FALSE)
        } else {
          q <- car::qqPlot(tpr, distribution = "chisq", df = ndim, 
                 lwd = 1, grid = FALSE, col="red",col.lines="black", main = "Multi-normal Q-Q Plot", xlab = expression(chi^2 * " quantiles"), ylab = expression("Mahalanobis distances "^2))
        }
      } else {
        qq <- lapply(pr,function(x) stats::qqnorm(x[[paste0(i,"Res")]],plot=FALSE))
        limInf <- min(min(unlist(lapply(qq,function(x) x$x)),na.rm=TRUE),min(unlist(lapply(qq,function(x) x$y)),na.rm=TRUE))
        limSup <- max(max(unlist(lapply(qq,function(x) x$x)),na.rm=TRUE),max(unlist(lapply(qq,function(x) x$y)),na.rm=TRUE))
        if (!requireNamespace("car", quietly = TRUE)) {
          q <- stats::qqnorm(pr[[1]][[paste0(i,"Res")]],main="",col="red",xlim=c(limInf,limSup),ylim=c(limInf,limSup))
          warning("Package \"car\" not installed. If pseudo-residual QQ plots using car::qqPlot are desired, please install it.",
               call. = FALSE)
        } else {
          q <- car::qqPlot(pr[[1]][[paste0(i,"Res")]],col="red",col.lines="black",main="",lwd=1,grid=FALSE,ylim=c(limInf,limSup),xlab="Theoretical Quantiles",ylab="Sample Quantiles")
        }
        if(nSims>1){
          for(j in 2:nSims){
            points(qq[[j]]$x,qq[[j]]$y,col=col[j])
          }
        }
        
        # add segments for zero inflation
        for(j in 1:nSims){
          if(m[[j]]$conditions$zeroInflation[[i]]) {
            ind <- which(m[[j]]$data[[i]]==0)
            x <- qq[[j]]$x[ind]
            y <- qq[[j]]$y[ind]
            segments(x,rep(limInf-5,length(ind)),x,y,col=col[j])
          }
          
          # add segments for one inflation
          if(m[[j]]$conditions$oneInflation[[i]]) {
            ind <- which(m[[j]]$data[[i]]==1)
            x <- qq[[j]]$x[ind]
            y <- qq[[j]]$y[ind]
            segments(x,rep(limInf-5,length(ind)),x,y,col=col[j])
          }
        }
        if (!requireNamespace("car", quietly = TRUE)) {
          abline(0,1,lwd=2)
        }
      }
      
      
      # ACF functions
      ACF <- lapply(pr,function(x) acf(x[[paste0(i,"Res")]],lag.max=lag.max,na.action=na.pass,plot=FALSE))
      acf(pr[[1]][[paste0(i,"Res")]],lag.max=lag.max,na.action=na.pass,main="",ylim=range(unlist(lapply(ACF,function(x) x$acf))))
      if(nSims>1){
        for(j in 2:nSims){
          lines(ACF[[j]]$lag,ACF[[j]]$acf,type="h",col=col[j])
        }
      }

  }

  # back to default
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # bottom, left, top, right
}

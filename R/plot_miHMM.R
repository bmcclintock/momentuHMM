
#' Plot \code{miHMM}
#'
#' Plot the fitted step and angle densities over histograms of the data, transition probabilities
#' as functions of the covariates, and maps of the animals' tracks colored by the decoded states.
#'
#' @method plot miHMM
#'
#' @param x Object \code{miHMM}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL} ; all animals are plotted.
#' @param covs Data frame consisting of a single row indicating the covariate values to be used in plots. If none are specified, the means of any covariates appearing in the model are used (unless covariate is a factor, in which case the first factor is used).
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param hist.ylim Parameter \code{ylim} for the step length histograms.
#' See \code{hist} documentation. Default: \code{NULL} ; the function sets default values.
#' @param sepAnimals If \code{TRUE}, the data is split by individuals in the histograms.
#' Default: \code{FALSE}.
#' @param sepStates If \code{TRUE}, the data is split by states in the histograms.
#' Default: \code{FALSE}.
#' @param col Vector or colors for the states (one color per state).
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @details The state-dependent densities are weighted by the frequency of each state in the most
#' probable state sequence (decoded with the function \code{\link{viterbi}}). For example, if the
#' most probable state sequence indicates that one third of observations correspond to the first
#' state, and two thirds to the second state, the plots of the densities in the first state are
#' weighted by a factor 1/3, and in the second state by a factor 2/3.
#'
#' @examples
#' # m is a miHMM object (as returned by MIfitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plot(m,ask=TRUE,animals=1,breaks=20)
#'
#'
#' @export
#' @importFrom graphics legend lines segments 
#' @importFrom grDevices adjustcolor

plot.miHMM <- function(x,animals=NULL,covs=NULL,ask=TRUE,breaks="Sturges",hist.ylim=NULL,sepAnimals=FALSE,
                              sepStates=FALSE,col=NULL,plotEllipse=TRUE,...)
{
  m <- x$miSum # the name "x" is for compatibility with the generic method
  nbAnimals <- length(unique(m$data$ID))
  nbStates <- length(m$stateNames)
  
  distnames <- names(m$conditions$dist)
  
  Fun <- lapply(m$conditions$dist,function(x) paste("d",x,sep=""))
  
  if(!is.null(hist.ylim) & length(hist.ylim)!=2)
    stop("hist.ylim needs to be a vector of two values (ymin,ymax)")
  
  # prepare colors for the states (used in the maps and for the densities)
  if(!is.null(col) & length(col)!=nbStates) {
    warning("Length of 'col' should be equal to number of states - argument ignored")
    col <- 2:(nbStates+1)
  }
  if(is.null(col))
    col <- 2:(nbStates+1)
  
  ######################
  ## Prepare the data ##
  ######################
  # determine indices of animals to be plotted
  if(is.null(animals)) # all animals are plotted
    animalsInd <- 1:nbAnimals
  else {
    if(is.character(animals)) { # animals' IDs provided
      animalsInd <- NULL
      for(zoo in 1:length(animals)) {
        if(length(which(unique(m$data$ID)==animals[zoo]))==0) # ID not found
          stop("Check 'animals' argument, ID not found")
        
        animalsInd <- c(animalsInd,which(unique(m$data$ID)==animals[zoo]))
      }
    }
    
    if(is.numeric(animals)) { # animals' indices provided
      if(length(which(animals<1))>0 | length(which(animals>nbAnimals))>0) # index out of bounds
        stop("Check 'animals' argument, index out of bounds")
      
      animalsInd <- animals
    }
  }
  
  nbAnimals <- length(animalsInd)
  ID <- unique(m$data$ID)[animalsInd]
  
  ##################################
  ## States decoding with Viterbi ##
  ##################################
  if(nbStates>1) {
    #cat("Decoding states sequence... ")
    states <- m$Par$states
    #cat("DONE\n")
  } else
    states <- rep(1,nrow(m$data))
  
  if(sepStates | nbStates==1)
    w <- rep(1,nbStates)
  else {
    # proportion of each state in the states sequence returned by the Viterbi algorithm
    w <- rep(NA,nbStates)
    for(state in 1:nbStates)
      w[state] <- length(which(states==state))/length(states)
  }
  
  if(all(c("x","y") %in% names(m$data))){
    x <- list()
    y <- list()
    errorEllipse <- list()
    for(zoo in 1:nbAnimals) {
      ind <- which(m$data$ID==ID[zoo])
      x[[zoo]] <- m$data$x[ind]
      y[[zoo]] <- m$data$y[ind]
      errorEllipse[[zoo]] <- m$errorEllipse[ind]
    }
  }
  
  #DMind <- lapply(m$conditions$fullDM,function(x) all(unlist(apply(x,1,function(y) lapply(y,length)))==1))
  #p <- parDef(m$conditions$dist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$DM,m$conditions$userBounds)
  if(is.null(covs)){
    covs <- m$data[which(m$data$ID %in% ID),][1,]
    for(j in names(m$data)[which(unlist(lapply(m$data,class))!="factor")]){
      covs[[j]]<-mean(m$data[[j]][which(m$data$ID %in% ID)],na.rm=TRUE)
    }
    #if(!is.null(m$rawCovs)){
    #  tmp <- apply(m$rawCovs,2,mean,na.rm=TRUE)
    #  covs <- matrix(tmp,nrow=1,ncol=length(tmp),dimnames=list("covs",names(tmp)))
    #  covs <- as.data.frame(covs)
    #} else covs <- m$data[1,]
  } else {
    if(!is.data.frame(covs)) stop('covs must be a data frame')
    if(nrow(covs)>1) stop('covs must consist of a single row')
    if(!all(names(covs) %in% names(m$data))) stop('invalid covs specified')
    if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(m$data$ID))
  }
  nbCovs <- ncol(model.matrix(m$conditions$formula,m$data))-1 # substract intercept column
  #par <- w2n(m$mod$estimate,m$conditions$bounds,p$parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,DMind,m$conditions$workcons,nrow(m$data),m$conditions$dist,p$Bndind)
  
  Par <- lapply(m$Par$real,function(x) x$est)
  #parindex <- c(0,cumsum(unlist(lapply(m$conditions$fullDM,ncol)))[-length(m$conditions$fullDM)])
  #names(parindex) <- distnames
  for(i in distnames){
    if(!is.null(m$conditions$DM[[i]])){# & m$conditions$DMind[[i]]){
      Par[[i]] <- m$Par$beta[[i]]$est#m$mod$estimate[parindex[[i]]+1:ncol(m$conditions$fullDM[[i]])]
      #names(Par[[i]])<-colnames(m$conditions$fullDM[[i]])
    }
  }
  Par<-lapply(Par,function(x) c(t(x)))
  Par<-Par[distnames]
  beta <- m$Par$beta$beta$est
  delta <- m$Par$real$delta$est
  
  for(i in distnames[which(m$conditions$dist %in% angledists)]){
    if(!m$conditions$estAngleMean[[i]]){
      m$conditions$estAngleMean[[i]]<-TRUE
      m$conditions$userBounds[[i]]<-rbind(matrix(rep(c(-pi,pi),nbStates),nbStates,2,byrow=TRUE),m$conditions$bounds[[i]])
      m$conditions$cons[[i]] <- c(rep(1,nbStates),m$conditions$cons[[i]])
      m$conditions$workcons[[i]] <- c(rep(0,nbStates),m$conditions$workcons[[i]])
      if(!is.null(m$conditions$DM[[i]])){
        Par[[i]] <- c(rep(0,nbStates),Par[[i]])
        if(is.list(m$conditions$DM[[i]])){
          m$conditions$DM[[i]]$mean<- ~1
        } else {
          tmpDM <- matrix(0,nrow(m$conditions$DM[[i]])+nbStates,ncol(m$conditions$DM[[i]])+nbStates)
          tmpDM[nbStates+1:nrow(m$conditions$DM[[i]]),nbStates+1:ncol(m$conditions$DM[[i]])] <- m$conditions$DM[[i]]
          diag(tmpDM)[1:nbStates] <- 1
          m$conditions$DM[[i]] <- tmpDM
        }
      }
    }
  }
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  p <- inputs$p
  DMinputs<-getDM(covs,inputs$DM,m$conditions$dist,nbStates,p$parNames,p$bounds,Par,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation)
  fullDM <- DMinputs$fullDM
  DMind <- DMinputs$DMind
  wpar <- n2w(Par,p$bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  par <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,inputs$estAngleMean,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,1,m$conditions$dist,p$Bndind)
  
  zeroMass<-vector('list',length(m$conditions$dist))
  names(zeroMass)<-distnames
  
  # text for legends
  stateNames <- m$stateNames
  if(is.null(stateNames)){
    legText <- NULL
    for(i in 1:nbStates)
      legText <- c(legText,paste("State",i))
  } else legText <- stateNames
  
  for(i in distnames){
    
    # split data by animals if necessary
    if(sepAnimals) {
      genData <- list()
      for(zoo in 1:nbAnimals) {
        ind <- which(m$data$ID==ID[zoo])
        genData[[zoo]] <- m$data[[i]][ind]
      }
    } else {
      ind <- which(m$data$ID %in% ID)
      genData <- m$data[[i]][ind]
    }
    
    if(m$conditions$zeroInflation[[i]]) {
      zeroMass[[i]] <- par[[i]][nrow(par[[i]])-(nbStates-1):0,]#m$mle[[i]][nrow(m$mle[[i]]),]
      par[[i]] <- par[[i]][-(nrow(par[[i]])-(nbStates-1):0),,drop=FALSE]#m$mle[[i]][-nrow(m$mle[[i]]),]
    } else {
      zeroMass[[i]] <- rep(0,nbStates)
    }
    
    infInd <- FALSE
    if(m$conditions$dist[[i]] %in% angledists)
      if(i=="angle" & ("step" %in% distnames))
        if(m$conditions$dist$step %in% stepdists & m$conditions$zeroInflation$step)
          infInd <- TRUE
    
    ###########################################
    ## Compute estimated densities on a grid ##
    ###########################################
    genDensities <- list()
    genFun <- Fun[[i]]
    if(m$conditions$dist[[i]] %in% angledists) {
      grid <- seq(-pi,pi,length=1000)
    } else if(m$conditions$dist[[i]]=="pois"){
      grid <- seq(0,max(m$data[[i]],na.rm=TRUE))
    } else {
      grid <- seq(0,max(m$data[[i]],na.rm=TRUE),length=10000)
    }
    
    for(state in 1:nbStates) {
      genArgs <- list(grid)
      
      for(j in 1:(nrow(par[[i]])/nbStates))
        genArgs[[j+1]] <- par[[i]][(j-1)*nbStates+state,]
      
      # conversion between mean/sd and shape/scale if necessary
      if(m$conditions$dist[[i]]=="gamma") {
        shape <- genArgs[[2]]^2/genArgs[[3]]^2
        scale <- genArgs[[3]]^2/genArgs[[2]]
        genArgs[[2]] <- shape
        genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
      }
      # (weighted by the proportion of each state in the Viterbi states sequence)
      if(m$conditions$zeroInflation[[i]]){
        genDensities[[state]] <- cbind(grid,(1-zeroMass[[i]][state])*w[state]*do.call(genFun,genArgs))
      } else if(infInd) {
        genDensities[[state]] <- cbind(grid,(1-zeroMass$step[state])*w[state]*do.call(genFun,genArgs))
      } else {
        genDensities[[state]] <- cbind(grid,w[state]*do.call(genFun,genArgs))
      }
    }
    
    #########################
    ## Plot the histograms ##
    #########################
    # set graphical parameters
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
    par(ask=ask)
    
    if(sepAnimals) {
      
      # loop over the animals
      for(zoo in 1:nbAnimals) {
        if(sepStates) {
          
          # loop over the states
          for(state in 1:nbStates) {
            gen <- genData[[zoo]][which(states[which(m$data$ID==ID[zoo])]==state)]
            message <- paste("Animal ID:",ID[zoo]," - State:",state)
            
            # the function plotHist is defined below
            plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,state,hist.ylim,col,legText)
          }
          
        } else { # if !sepStates
          gen <- genData[[zoo]]
          message <- paste("Animal ID:",ID[zoo])
          
          plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,NULL,hist.ylim,col,legText)
        }
      }
    } else { # if !sepAnimals
      if(sepStates) {
        
        # loop over the states
        for(state in 1:nbStates) {
          gen <- genData[which(states==state)]
          message <- paste("All animals - State:",state)
          
          plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,state,hist.ylim,col,legText)
        }
        
      } else { # if !sepStates
        gen <- genData
        message <- "All animals"
        
        plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,NULL,hist.ylim,col,legText)
      }
    }
  }
  
  ##################################################
  ## Plot the t.p. as functions of the covariates ##
  ##################################################
  if(nbStates>1) {
    par(mfrow=c(nbStates,nbStates))
    par(mar=c(5,4,4,2)-c(0,0,1.5,1)) # bottom, left, top, right
    
    rawCovs <- m$rawCovs
    gridLength <- 100
    
    if(nrow(beta)>1) {
      for(cov in 1:ncol(m$rawCovs)) {
        inf <- min(rawCovs[,cov],na.rm=T)
        sup <- max(rawCovs[,cov],na.rm=T)
        
        # mean values of each covariate
        meanCovs <- colSums(rawCovs)/nrow(rawCovs)
        
        # set all covariates to their mean, except for "cov"
        # (which takes a grid of values from inf to sup)
        tempCovs <- data.frame(rep(meanCovs[1],gridLength))
        if(length(meanCovs)>1)
          for(i in 2:length(meanCovs))
            tempCovs <- cbind(tempCovs,rep(meanCovs[i],gridLength))
        
        tempCovs[,cov] <- seq(inf,sup,length=gridLength)
        colnames(tempCovs) <- colnames(rawCovs)
        
        desMat <- model.matrix(m$conditions$formula,data=tempCovs)
        
        # check that the current covariate (cov) is included in the model
        used <- FALSE
        for(i in 2:ncol(desMat)) {
          c <- desMat[,i]
          if(length(which(c!=mean(c)))>0)
            used <- TRUE
        }
        
        if(used) {
          trMat <- trMatrix_rcpp(nbStates,beta,desMat)
          
          for(i in 1:nbStates)
            for(j in 1:nbStates)
              plot(tempCovs[,cov],trMat[i,j,],type="l",ylim=c(0,1),xlab=names(rawCovs)[cov],
                   ylab=paste(i,"->",j))
          
          mtext("Transition probabilities",side=3,outer=TRUE,padj=2)
        }
      }
    }
  }
  
  #################################
  ## Plot maps colored by states ##
  #################################
  if(all(c("x","y") %in% names(m$data))){
    
    if(nbStates>1) { # no need to plot the map if only one state
      par(mfrow=c(1,1))
      par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
      
      for(zoo in 1:nbAnimals) {
        # states for animal 'zoo'
        subStates <- states[which(m$data$ID==ID[zoo])]
        
        # determine the bounds of the plot
        xmin <- min(x[[zoo]],na.rm=T)
        xmax <- max(x[[zoo]],na.rm=T)
        ymin <- min(y[[zoo]],na.rm=T)
        ymax <- max(y[[zoo]],na.rm=T)
        # make sure that x and y have same scale
        if(xmax-xmin>ymax-ymin) {
          ymid <- (ymax+ymin)/2
          ymax <- ymid+(xmax-xmin)/2
          ymin <- ymid-(xmax-xmin)/2
        } else {
          xmid <- (xmax+xmin)/2
          xmax <- xmid+(ymax-ymin)/2
          xmin <- xmid-(ymax-ymin)/2
        }
        
        # first point
        plot(x[[zoo]][1],y[[zoo]][1],xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=18,
             xlab="x",ylab="y",col=col[subStates[1]])
        if(plotEllipse) lines(errorEllipse[[zoo]][[1]],col=adjustcolor(col[subStates[1]],alpha.f=0.25),cex=0.6)
        
        # trajectory
        for(i in 2:length(x[[zoo]])) {
          points(x[[zoo]][i],y[[zoo]][i],pch=16,col=col[subStates[i-1]],cex=0.6)
          segments(x0=x[[zoo]][i-1],y0=y[[zoo]][i-1],x1=x[[zoo]][i],y1=y[[zoo]][i],
                   col=col[subStates[i-1]],lwd=1.3)
          if(plotEllipse) lines(errorEllipse[[zoo]][[i]],col=adjustcolor(col[subStates[i-1]],alpha.f=0.25),cex=0.6)
        }
        mtext(paste("Animal ID:",ID[zoo]),side=3,outer=TRUE,padj=2)
        legend("topleft",legText,lwd=rep(2,nbStates),col=col,bty="n")
      }
    }
  }
  
  # set the graphical parameters back to default
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # bottom, left, top, right
  par(ask=FALSE)
}

# Plot histograms
#
# Plot histograms of steps and angles, and the fitted densities. This function is only
# used in the function plot.miHMM.
#
# Parameters:
#  - step: list of series of steps if several animals, or series of steps otherwise.
#    (e.g. step[[1]][3] is the 3rd step of the first animal)
#  - angle: same as step, but for angles
#  - stepDensities: list of matrices of values of the fitted densities. Each matrix has
#    two columns, the first being the grid of values on which the density is estimated,
#    and the second the values of the density.
#  - angleDensities: same as stepDensities, but for angles.
#  - message: message to print above the histograms
#  - sepStates, breaks, hist.ylim: see arguments of plot.miHMM.
#  - state: if sepStates, this function needs to know which state needs to be plotted.
#  - col: colors of the state-dependent density lines

plotHist <- function (gen,genDensities,dist,message,
                      sepStates,breaks="Sturges",state=NULL,hist.ylim=NULL,col=NULL,legText)
{
  # vertical limits
  if(!is.null(hist.ylim)) {
    ymin <- hist.ylim[1]
    ymax <- hist.ylim[2]
  } else {
    ymin <- 0
    ymax <- NA
  }
  
  if(!sepStates) {
    nbStates <- length(genDensities)
  }
  
  distname <- names(dist)
  
  if(dist %in% angledists){
    h <- hist(gen,plot=F,breaks=breaks) # to determine 'breaks'
    breaks <- seq(-pi,pi,length=length(h$breaks))
    
    h <- hist(gen,plot=F,breaks=breaks) # to determine 'ymax'
    ymax <- 1.3*max(h$density)
    
    # plot gen histogram
    hist(gen,prob=T,main="",ylim=c(0,ymax),xlab=paste0(distname," (radians)"),
         col="grey",border="white",breaks=breaks,xaxt="n")
    axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
         labels = expression(-pi, -pi/2, 0, pi/2, pi))
    
    mtext(message,side=3,outer=TRUE,padj=2)
    
    # plot gen density over the histogram
    if(sepStates)
      lines(genDensities[[state]],col=col[state],lwd=2)
    else {
      for(s in 1:nbStates)
        lines(genDensities[[s]],col=col[s],lwd=2)
      
      legend("top",legText,lwd=rep(2,nbStates),col=col,bty="n")
    }  
  } else {
    # determine ylim
    if(is.null(hist.ylim)) { # default
      h <- hist(gen,plot=F,breaks=breaks)
      ymax <- 1.3*max(h$density)
      
      # find the maximum of the gen densit-y-ies, and take it as ymax if necessary
      if(sepStates) {
        maxdens <- max(genDensities[[state]][,2])
        if(maxdens>ymax & maxdens<2*max(h$density))
          ymax <- maxdens
        
      } else {
        maxdens <- max(genDensities[[1]][,2])
        if(nbStates>1) {
          for(state in 2:nbStates) {
            if(is.finite(max(genDensities[[state]][,2]))){
              if(max(genDensities[[state]][,2])>maxdens)
                maxdens <- max(genDensities[[state]][,2])
            }
          }
        }
        if(maxdens>ymax & maxdens<2*max(h$density))
          ymax <- maxdens
      }
    }
    
    # plot gen histogram
    hist(gen,prob=T,main="",ylim=c(ymin,ymax),xlab=distname,
         col="grey",border="white",breaks=breaks)
    
    mtext(message,side=3,outer=TRUE,padj=2)
    
    # plot gen density over the histogram
    if(sepStates)
      lines(genDensities[[state]],col=col[state],lwd=2)
    else {
      for(s in 1:nbStates)
        lines(genDensities[[s]],col=col[s],lwd=2)
      
      legend("top",legText,lwd=rep(2,nbStates),col=col,bty="n")
    }
  }
  
}

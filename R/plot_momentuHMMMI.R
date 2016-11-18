
#' Plot \code{momentuHMMMI}
#'
#' Plot the fitted step and angle densities over histograms of the data, transition probabilities
#' as functions of the covariates, and maps of the animals' tracks colored by the decoded states.
#'
#' @method plot momentuHMMMI
#'
#' @param x Object \code{momentuHMMMI}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL} ; all animals are plotted.
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
#'
#' @export
#' @importFrom graphics legend lines segments

plot.momentuHMMMI <- function(x,animals=NULL,ask=TRUE,breaks="Sturges",hist.ylim=NULL,sepAnimals=FALSE,
                         sepStates=FALSE,col=NULL,...)
{
  m <- x # the name "x" is for compatibility with the generic method
  nbAnimals <- length(unique(m$data$ID))
  nbStates <- ncol(m$Par$stepPar$est)

  stepFun <- paste("d",m$conditions$stepDist,sep="")
  if(m$conditions$angleDist!="none")
    angleFun <- paste("d",m$conditions$angleDist,sep="")

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

  # split data by animals if necessary
  if(sepAnimals) {
    stepData <- list()
    angleData <- list()
    for(zoo in 1:nbAnimals) {
      ind <- which(m$data$ID==ID[zoo])
      stepData[[zoo]] <- m$data$step[ind]
      angleData[[zoo]] <- m$data$angle[ind]
    }
  } else {
    ind <- which(m$data$ID %in% ID)
    stepData <- m$data$step[ind]
    angleData <- m$data$angle[ind]
  }

  if(m$conditions$angleDist=="none")
    angleData <- NULL

  x <- list()
  y <- list()
  for(zoo in 1:nbAnimals) {
    ind <- which(m$data$ID==ID[zoo])
    x[[zoo]] <- m$data$x[ind]
    y[[zoo]] <- m$data$y[ind]
  }

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

  if(m$conditions$zeroInflation) {
    zeromass <- m$Par$stepPar$est[nrow(m$Par$stepPar$est),]
    m$Par$stepPar <- m$Par$stepPar$est[-nrow(m$Par$stepPar$est),]
  }

  ###########################################
  ## Compute estimated densities on a grid ##
  ###########################################
  stepDensities <- list()
  grid <- seq(0,max(m$data$step,na.rm=TRUE),length=10000)

  for(state in 1:nbStates) {
    stepArgs <- list(grid)

    for(j in 1:nrow(m$Par$stepPar$est))
      stepArgs[[j+1]] <- m$Par$stepPar$est[j,state]

    # conversion between mean/sd and shape/scale if necessary
    if(m$conditions$stepDist=="gamma") {
      shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
      scale <- stepArgs[[3]]^2/stepArgs[[2]]
      stepArgs[[2]] <- shape
      stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
    }
    # (weighted by the proportion of each state in the Viterbi states sequence)
    if(m$conditions$zeroInflation)
      stepDensities[[state]] <- cbind(grid,(1-zeromass[state])*w[state]*do.call(stepFun,stepArgs))
    else
      stepDensities[[state]] <- cbind(grid,w[state]*do.call(stepFun,stepArgs))
  }

  if(m$conditions$angleDist!="none") {
    angleDensities <- list()
    grid <- seq(-pi,pi,length=1000)

    for(state in 1:nbStates) {
      angleArgs <- list(grid)

      for(j in 1:nrow(m$Par$anglePar$est))
        angleArgs[[j+1]] <- m$Par$anglePar$est[j,state]

      # (weighted by the proportion of each state in the Viterbi states sequence)
      angleDensities[[state]] <- cbind(grid,w[state]*do.call(angleFun,angleArgs))
    }
  }
  
  # text for legends
  stateNames <- m$stateNames
  if(is.null(stateNames)){
    legText <- NULL
    for(i in 1:nbStates)
      legText <- c(legText,paste("State",i))
  } else legText <- stateNames

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
          step <- stepData[[zoo]][which(states[which(m$data$ID==ID[zoo])]==state)]
          angle <- angleData[[zoo]][which(states[which(m$data$ID==ID[zoo])]==state)]
          message <- paste("Animal ID:",ID[zoo]," - State:",state)

          # the function plotHist is defined below
          plotHist(step,angle,stepDensities,angleDensities,message,sepStates,breaks,state,hist.ylim,col,legText)
        }

      } else { # if !sepStates
        step <- stepData[[zoo]]
        angle <- angleData[[zoo]]
        message <- paste("Animal ID:",ID[zoo])

        plotHist(step,angle,stepDensities,angleDensities,message,sepStates,breaks,NULL,hist.ylim,col,legText)
      }
    }
  } else { # if !sepAnimals
    if(sepStates) {

      # loop over the states
      for(state in 1:nbStates) {
        step <- stepData[which(states==state)]
        angle <- angleData[which(states==state)]
        message <- paste("All animals - State:",state)

        plotHist(step,angle,stepDensities,angleDensities,message,sepStates,breaks,state,hist.ylim,col,legText)
      }

    } else { # if !sepStates
      step <- stepData
      angle <- angleData
      message <- "All animals"

      plotHist(step,angle,stepDensities,angleDensities,message,sepStates,breaks,NULL,hist.ylim,col,legText)
    }
  }

  ##################################################
  ## Plot the t.p. as functions of the covariates ##
  ##################################################
  if(!is.null(m$Par$beta$est) & nbStates>1) {
    par(mfrow=c(nbStates,nbStates))
    par(mar=c(5,4,4,2)-c(0,0,1.5,1)) # bottom, left, top, right

    rawCovs <- m$rawCovs
    gridLength <- 100

    if(nrow(m$Par$beta$est)>1) {
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
          trMat <- trMatrix_rcpp(nbStates,m$Par$beta$est,desMat)

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

      # trajectory
      for(i in 2:length(x[[zoo]])) {
        points(x[[zoo]][i],y[[zoo]][i],pch=16,col=col[subStates[i-1]],cex=0.6)
        segments(x0=x[[zoo]][i-1],y0=y[[zoo]][i-1],x1=x[[zoo]][i],y1=y[[zoo]][i],
                 col=col[subStates[i-1]],lwd=1.3)
      }
      mtext(paste("Animal ID:",ID[zoo]),side=3,outer=TRUE,padj=2)
      legend("topleft",legText,lwd=rep(2,nbStates),col=col,bty="n")
    }
  }

  # set the graphical parameters back to default
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # bottom, left, top, right
  par(ask=FALSE)
}
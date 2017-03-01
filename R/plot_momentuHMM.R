
#' Plot \code{momentuHMM}
#'
#' Plot the fitted step and angle densities over histograms of the data, transition probabilities
#' as functions of the covariates, and maps of the animals' tracks colored by the decoded states.
#'
#' @method plot momentuHMM
#'
#' @param x Object \code{momentuHMM}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL} ; all animals are plotted.
#' @param covs Data frame consisting of a single row indicating the covariate values to be used in plots. 
#' If none are specified, the means of any covariates appearing in the model are used (unless covariate is a factor, in which case the first factor in the data is used).
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param hist.ylim An optional named list of vectors specifying \code{ylim=c(ymin,ymax)} for the data stream histograms.
#' See \code{hist} documentation. Default: \code{NULL} ; the function sets default values for all data streams.
#' @param sepAnimals If \code{TRUE}, the data is split by individuals in the histograms.
#' Default: \code{FALSE}.
#' @param sepStates If \code{TRUE}, the data is split by states in the histograms.
#' Default: \code{FALSE}.
#' @param col Vector or colors for the states (one color per state).
#' @param cumul	If TRUE, the sum of weighted densities is plotted (default).
#' @param plotTracks If TRUE, the Viterbi-decoded tracks are plotted (default).
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @details The state-dependent densities are weighted by the frequency of each state in the most
#' probable state sequence (decoded with the function \code{\link{viterbi}}). For example, if the
#' most probable state sequence indicates that one third of observations correspond to the first
#' state, and two thirds to the second state, the plots of the densities in the first state are
#' weighted by a factor 1/3, and in the second state by a factor 2/3.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plot(m,ask=TRUE,animals=1,breaks=20)
#'
#'
#' @export
#' @importFrom graphics legend lines segments arrows
#' @importFrom grDevices adjustcolor gray rainbow
#' @importFrom stats as.formula

plot.momentuHMM <- function(x,animals=NULL,covs=NULL,ask=TRUE,breaks="Sturges",hist.ylim=NULL,sepAnimals=FALSE,
                         sepStates=FALSE,col=NULL,cumul=TRUE,plotTracks=TRUE,plotCI=FALSE,alpha=0.95,...)
{
  m <- x # the name "x" is for compatibility with the generic method
  nbAnimals <- length(unique(m$data$ID))
  nbStates <- length(m$stateNames)
  
  distnames <- names(m$conditions$dist)

  Fun <- lapply(m$conditions$dist,function(x) paste("d",x,sep=""))

  if(is.null(hist.ylim)){
    hist.ylim<-vector('list',length(distnames))
    names(hist.ylim)<-distnames
  }
  for(i in distnames){
    if(!is.null(hist.ylim[[i]]) & length(hist.ylim[[i]])!=2)
      stop("hist.ylim$",i," needs to be a vector of two values (ymin,ymax)")
  }
  
  # prepare colors for the states (used in the maps and for the densities)
  if (!is.null(col) & length(col) != nbStates) {
    warning("Length of 'col' should be equal to number of states - argument ignored")
    col <- 2:(nbStates + 1)
  }
  if (is.null(col) & nbStates < 8) {
    pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
             "#0072B2", "#D55E00", "#CC79A7")
    col <- pal[1:nbStates]
  }
  if (is.null(col) & nbStates >= 8) 
    col <- rainbow(nbStates)
  
  if (sepStates | nbStates < 2) 
    cumul <- FALSE
  
  if("miSum" %in% class(x)) plotEllipse <- m$plotEllipse
  else plotEllipse <- FALSE
  
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
    if("miSum" %in% class(x)) states <- m$Par$states
    else {
      cat("Decoding states sequence... ")
      states <- viterbi(m)
      cat("DONE\n")
    }
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
    if(plotEllipse)  errorEllipse <- list()
    for(zoo in 1:nbAnimals) {
      ind <- which(m$data$ID==ID[zoo])
      x[[zoo]] <- m$data$x[ind]
      y[[zoo]] <- m$data$y[ind]
      if(plotEllipse) errorEllipse[[zoo]] <- m$errorEllipse[ind]
    }
  }
  
  if(is.null(covs)){
    covs <- m$data[which(m$data$ID %in% ID),][1,]
    for(j in names(m$data)[which(unlist(lapply(m$data,class))!="factor")]){
      covs[[j]]<-mean(m$data[[j]][which(m$data$ID %in% ID)],na.rm=TRUE)
    }
  } else {
    if(!is.data.frame(covs)) stop('covs must be a data frame')
    if(nrow(covs)>1) stop('covs must consist of a single row')
    if(!all(names(covs) %in% names(m$data))) stop('invalid covs specified')
    if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(m$data$ID))
    for(j in names(m$data)[which(!(names(m$data) %in% names(covs)))]){
      if(class(m$data[[j]])=="factor") covs[[j]] <- m$data[[j]][which(m$data$ID %in% ID)][1]
      else covs[[j]]<-mean(m$data[[j]][which(m$data$ID %in% ID)],na.rm=TRUE)
    }
    for(j in names(m$data)[which(names(m$data) %in% names(covs))]){
      if(class(m$data[[j]])=="factor") covs[[j]] <- factor(covs[[j]],levels=levels(m$data[[j]]))
      if(is.na(covs[[j]])) stop("check value for ",j)
    }
  }
  nbCovs <- ncol(model.matrix(m$conditions$formula,m$data))-1 # substract intercept column
 
  Par <- m$mle[distnames]
  parindex <- c(0,cumsum(unlist(lapply(m$conditions$fullDM,ncol)))[-length(m$conditions$fullDM)])
  names(parindex) <- distnames
  for(i in distnames){
    if(!is.null(m$conditions$DM[[i]])){# & m$conditions$DMind[[i]]){
      Par[[i]] <- m$mod$estimate[parindex[[i]]+1:ncol(m$conditions$fullDM[[i]])]
      names(Par[[i]])<-colnames(m$conditions$fullDM[[i]])
    }
  }
  Par<-lapply(Par,function(x) c(t(x)))
  beta <- m$mle$beta
  delta <- m$mle$delta
  
  tmpPar <- Par
  tmpConditions <- m$conditions
  
  for(i in distnames[which(m$conditions$dist %in% angledists)]){
    if(!m$conditions$estAngleMean[[i]]){
      tmpConditions$estAngleMean[[i]]<-TRUE
      tmpConditions$userBounds[[i]]<-rbind(matrix(rep(c(-pi,pi),nbStates),nbStates,2,byrow=TRUE),m$conditions$bounds[[i]])
      tmpConditions$cons[[i]] <- c(rep(1,nbStates),m$conditions$cons[[i]])
      tmpConditions$workcons[[i]] <- c(rep(0,nbStates),m$conditions$workcons[[i]])
      if(!is.null(m$conditions$DM[[i]])){
        tmpPar[[i]] <- c(rep(0,nbStates),Par[[i]])
        if(is.list(m$conditions$DM[[i]])){
          tmpConditions$DM[[i]]$mean<- ~1
        } else {
          tmpDM <- matrix(0,nrow(tmpConditions$DM[[i]])+nbStates,ncol(tmpConditions$DM[[i]])+nbStates)
          tmpDM[nbStates+1:nrow(tmpConditions$DM[[i]]),nbStates+1:ncol(tmpConditions$DM[[i]])] <- tmpConditions$DM[[i]]
          diag(tmpDM)[1:nbStates] <- 1
          tmpConditions$DM[[i]] <- tmpDM
        }
      } else {
        Par[[i]] <- Par[[i]][-(1:nbStates)]
      }
    }
  }
  
  # get pars for probability density plots
  tmpInputs <- checkInputs(nbStates,tmpConditions$dist,tmpPar,tmpConditions$estAngleMean,tmpConditions$circularAngleMean,tmpConditions$zeroInflation,tmpConditions$DM,tmpConditions$userBounds,tmpConditions$cons,tmpConditions$workcons,m$stateNames)
  tmpp <- tmpInputs$p
  DMinputs<-getDM(covs,tmpInputs$DM,tmpConditions$dist,nbStates,tmpp$parNames,tmpp$bounds,tmpPar,tmpConditions$cons,tmpConditions$workcons,tmpConditions$zeroInflation,tmpConditions$circularAngleMean)
  fullDM <- DMinputs$fullDM
  DMind <- DMinputs$DMind
  wpar <- n2w(tmpPar,tmpp$bounds,beta,rep(1/nbStates,nbStates),nbStates,tmpInputs$estAngleMean,tmpInputs$DM,DMinputs$cons,DMinputs$workcons,tmpp$Bndind)
  if(!m$conditions$stationary & nbStates>1) {
    wpar[(length(wpar)-nbStates+2):length(wpar)] <- m$mod$estimate[(length(m$mod$estimate)-nbStates+2):length(m$mod$estimate)] #this is done to deal with any delta=0 in n2w
  }
  par <- w2n(wpar,tmpp$bounds,tmpp$parSize,nbStates,nbCovs,tmpInputs$estAngleMean,tmpInputs$circularAngleMean,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,1,tmpConditions$dist,tmpp$Bndind)
  
  inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  p <- inputs$p
  
  zeroMass<-vector('list',length(m$conditions$dist))
  names(zeroMass)<-distnames
  
  # text for legends
  stateNames <- m$stateNames
  if(is.null(stateNames)){
    legText <- NULL
    for(i in 1:nbStates)
      legText <- c(legText,paste("State",i))
  } else legText <- stateNames
  
  tmpcovs<-covs
  for(i in which(!mapply(is.factor,covs))){
    tmpcovs[i]<-round(covs[i],2)
  }
  for(i in which(mapply(is.factor,covs))){
    tmpcovs[i] <- as.character(covs[[i]])
  }
  
  Sigma <- ginv(m$mod$hessian)
  
  # set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
  par(ask=ask)
  
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
          if(all(c("x","y") %in% names(m$data)))
            infInd <- TRUE
    
    #get covariate names
    DMparterms<-list()
    if(!m$conditions$DMind[[i]]){
      if(!is.list(m$conditions$DM[[i]])){
        for(j in 1:length(p$parNames[[i]])){
          DMparterms[[p$parNames[[i]][j]]] <- vector('list',nbStates)
          for(jj in 1:nbStates){
            DMparterms[[p$parNames[[i]][j]]][[jj]]<-unique(m$conditions$DM[[i]][(j-1)*nbStates+jj,][suppressWarnings(which(is.na(as.numeric(m$conditions$DM[[i]][(j-1)*nbStates+jj,]))))])
          }
        }
        DMterms<-unique(m$conditions$DM[[i]][suppressWarnings(which(is.na(as.numeric(m$conditions$DM[[i]]))))])
      } else {
        m$conditions$DM[[i]]<-m$conditions$DM[[i]][p$parNames[[i]]]
        DMterms<-character()
        for(j in 1:length(p$parNames[[i]])){
          DMparterms[[p$parNames[[i]][j]]] <- vector('list',nbStates)
          tmpparnames <- rownames(attr(terms(m$conditions$DM[[i]][[p$parNames[[i]][j]]]),"factors"))
          if(!is.null(tmpparnames)) {
            for(jj in 1:nbStates){
              DMparterms[[p$parNames[[i]][j]]][[jj]]<-tmpparnames
            }
          }
          DMterms <- c(DMterms,unlist(DMparterms[[p$parNames[[i]][j]]]))
        }
      }
      DMterms <- unique(DMterms)
      for(j in 1:length(p$parNames[[i]])){
        for(jj in 1:nbStates){
          if(length(DMparterms[[p$parNames[[i]][j]]][[jj]])){
            for(k in 1:length(DMparterms[[p$parNames[[i]][j]]][[jj]])){
              DMparterms[[p$parNames[[i]][j]]][[jj]][k]<-all.vars(as.formula(paste0("~",DMparterms[[p$parNames[[i]][j]]][[jj]][k])))
            }
          }
        }
      }
      for(j in 1:length(DMterms)){
        DMterms[j]<-all.vars(as.formula(paste0("~",DMterms[j])))
      }
    }
    
    covmess <- ifelse(!m$conditions$DMind[[i]],paste0(": ",paste0(DMterms," = ",tmpcovs[DMterms],collapse=", ")),"")
  
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
      
      for(j in p$parNames[[i]]){
      
        for(jj in DMparterms[[j]][[state]]){
          
          if(!is.factor(m$data[,jj])){
            
            gridLength <- 101
            
            inf <- min(m$data[,jj],na.rm=T)
            sup <- max(m$data[,jj],na.rm=T)
            
            # set all covariates to their mean, except for "cov"
            # (which takes a grid of values from inf to sup)
            tempCovs <- data.frame(matrix(covs[DMparterms[[j]][[state]]][[1]],nrow=gridLength,ncol=1))
            if(length(DMterms)>1)
              for(ii in 2:length(DMterms))
                tempCovs <- cbind(tempCovs,rep(covs[[DMterms[[ii]]]],gridLength))
            names(tempCovs) <- DMterms
            tempCovs[,jj] <- seq(inf,sup,length=gridLength)
          } else {
            gridLength<- nlevels(m$data[,jj])
            # set all covariates to their mean, except for "cov"
            tempCovs <- data.frame(matrix(covs[DMparterms[[j]][[state]]][[1]],nrow=gridLength,ncol=1))
            if(length(DMterms)>1)
              for(ii in 2:length(DMterms))
                tempCovs <- cbind(tempCovs,rep(covs[[DMterms[[ii]]]],gridLength))
            names(tempCovs) <- DMterms
            tempCovs[,jj] <- as.factor(levels(m$data[,jj]))
          }
          
          for(ii in DMterms[which(unlist(lapply(m$data[DMterms],is.factor)))])
            tempCovs[[ii]] <- factor(tempCovs[[ii]],levels=levels(m$data[[ii]]))
          
          DMinputs<-getDM(tempCovs,inputs$DM[i],m$conditions$dist[i],nbStates,p$parNames[i],p$bounds[i],Par[i],m$conditions$cons[i],m$conditions$workcons[i],m$conditions$zeroInflation[i],m$conditions$circularAngleMean[i])
          fullDM <- DMinputs$fullDM
          DMind <- DMinputs$DMind
          gradfun<-function(wpar,k) {
            w2n(wpar,p$bounds[i],p$parSize[i],nbStates,nbCovs,inputs$estAngleMean[i],inputs$circularAngleMean[i],stationary=TRUE,DMinputs$cons[i],fullDM,DMind,DMinputs$workcons[i],gridLength,m$conditions$dist[i],p$Bndind[i])[[i]][(which(p$parNames[[i]]==j)-1)*nbStates+state,k]
          }
          est<-w2n(c(m$mod$estimate[parindex[[i]]+1:ncol(fullDM[[i]])],beta),p$bounds[i],p$parSize[i],nbStates,nbCovs,inputs$estAngleMean[i],inputs$circularAngleMean[i],stationary=TRUE,DMinputs$cons[i],fullDM,DMind,DMinputs$workcons[i],gridLength,m$conditions$dist[i],p$Bndind[i])[[i]][(which(p$parNames[[i]]==j)-1)*nbStates+state,]
          if(plotCI){
            dN<-t(mapply(function(x) tryCatch(numDeriv::grad(gradfun,c(m$mod$estimate[parindex[[i]]+1:ncol(fullDM[[i]])],beta),k=x),error=function(e) NA),1:gridLength))
            se<-t(apply(dN[,1:ncol(fullDM[[i]])],1,function(x) tryCatch(suppressWarnings(sqrt(x%*%Sigma[parindex[[i]]+1:ncol(fullDM[[i]]),parindex[[i]]+1:ncol(fullDM[[i]])]%*%x)),error=function(e) NA)))
            uci<-est+qnorm(1-(1-alpha)/2)*se
            lci<-est-qnorm(1-(1-alpha)/2)*se
            plot(tempCovs[,jj],est,ylim=range(c(lci,est,uci),na.rm=TRUE),xaxt="n",xlab=jj,ylab=paste(i,j,'parameter'),main=paste0('State ',state,ifelse(length(tempCovs[,DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]]==j)]]),paste0(": ",paste(DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]]==j)],"=",tmpcovs[,DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]]==j)]],collapse=", ")),"")),type="l")
            if(!all(is.na(se))){
              seInd <- which(!is.na(se))
              ciInd <- seInd[which(abs(uci[seInd]-lci[seInd])>max(abs(uci[seInd]-lci[seInd]))/1000)] #to aviod un-supressable warning in arrows()
              arrows(as.numeric(tempCovs[ciInd,jj]), lci[ciInd], as.numeric(tempCovs[ciInd,jj]), uci[ciInd], length=0.025, angle=90, code=3, col=gray(.5)) 
            }
          } else plot(tempCovs[,jj],est,xaxt="n",xlab=jj,ylab=paste(i,j,'parameter'),main=paste0('State ',state,ifelse(length(tempCovs[,DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]]==j)]]),paste0(": ",paste(DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]]==j)],"=",tmpcovs[,DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]]==j)]],collapse=", ")),"")),type="l") 
          if(is.factor(tempCovs[,jj])) axis(1,at=tempCovs[,jj])
          else axis(1)
          
        }
      }
      
    }
  
    #########################
    ## Plot the histograms ##
    #########################
    if(sepAnimals) {
  
      # loop over the animals
      for(zoo in 1:nbAnimals) {
        if(sepStates) {
  
          # loop over the states
          for(state in 1:nbStates) {
            gen <- genData[[zoo]][which(states[which(m$data$ID==ID[zoo])]==state)]
            message <- paste0("Animal ID ",ID[zoo]," - State ",state,covmess)
  
            # the function plotHist is defined below
            plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,state,hist.ylim[[i]],col,legText, cumul = cumul)
          }
  
        } else { # if !sepStates
          gen <- genData[[zoo]]
          message <- paste0("Animal ID ",ID[zoo],covmess)
  
          plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,NULL,hist.ylim[[i]],col,legText, cumul = cumul)
        }
      }
    } else { # if !sepAnimals
      if(sepStates) {
  
        # loop over the states
        for(state in 1:nbStates) {
          gen <- genData[which(states==state)]
          message <- paste0("All animals - State ",state,covmess)
  
          plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,state,hist.ylim[[i]],col,legText, cumul = cumul)
        }
  
      } else { # if !sepStates
        gen <- genData
        message <- paste0("All animals",covmess)
  
        plotHist(gen,genDensities,m$conditions$dist[i],message,sepStates,breaks,NULL,hist.ylim[[i]],col,legText, cumul = cumul)
      }
    }
  }

  ##################################################
  ## Plot the t.p. as functions of the covariates ##
  ##################################################
  if(nbStates>1) {
    par(mfrow=c(nbStates,nbStates))
    par(mar=c(5,4,4,2)-c(0,0,1.5,1)) # bottom, left, top, right
    
    gamInd<-(length(m$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)+1):(length(m$mod$estimate))-(nbStates-1)
    quantSup<-qnorm(1-(1-alpha)/2)
    
    if(nrow(m$mle$beta)>1) {
      
      # values of each covariate
      rawCovs <- m$rawCovs[which(m$data$ID %in% ID),,drop=FALSE]
      #if(is.null(covs)) {
      #  rawCovs <- m$rawCovs
      #  meanCovs <- colSums(rawCovs)/nrow(rawCovs)
      #} else {
      #  rawCovs <- m$data[,names(covs),drop=FALSE]
      #  meanCovs <- as.numeric(covs)
      #}
      
      for(cov in 1:ncol(rawCovs)) {
        
        if(!is.factor(rawCovs[,cov])){
          
          gridLength <- 101
          
          inf <- min(rawCovs[,cov],na.rm=T)
          sup <- max(rawCovs[,cov],na.rm=T)
          
          # set all covariates to their mean, except for "cov"
          # (which takes a grid of values from inf to sup)
          tempCovs <- data.frame(matrix(covs[names(rawCovs)][[1]],nrow=gridLength,ncol=1))
          if(ncol(rawCovs)>1)
            for(i in 2:ncol(rawCovs))
              tempCovs <- cbind(tempCovs,rep(covs[names(rawCovs)][[i]],gridLength))
          
          tempCovs[,cov] <- seq(inf,sup,length=gridLength)
        } else {
          gridLength<- nlevels(rawCovs[,cov])
          # set all covariates to their mean, except for "cov"
          tempCovs <- data.frame(matrix(covs[names(rawCovs)][[1]],nrow=gridLength,ncol=1))
          if(ncol(rawCovs)>1)
            for(i in 2:ncol(rawCovs))
              tempCovs <- cbind(tempCovs,rep(covs[names(rawCovs)][[i]],gridLength))
          
          tempCovs[,cov] <- as.factor(levels(rawCovs[,cov]))
        }
        
        names(tempCovs) <- names(rawCovs)
        tmpcovs<-covs[names(rawCovs)]
        for(i in which(unlist(lapply(rawCovs,is.factor)))){
          tempCovs[[i]] <- factor(tempCovs[[i]],levels=levels(rawCovs[,i]))
          tmpcovs[i] <- as.character(tmpcovs[[i]])
        }
        for(i in which(!unlist(lapply(rawCovs,is.factor)))){
          tmpcovs[i]<-round(covs[names(rawCovs)][i],2)
        }
        
        desMat <- model.matrix(m$conditions$formula,data=tempCovs)
        
        trMat <- trMatrix_rcpp(nbStates,m$mle$beta,desMat)
          
        for(i in 1:nbStates){
          for(j in 1:nbStates){
            plot(tempCovs[,cov],trMat[i,j,],type="l",ylim=c(0,1),xlab=names(rawCovs)[cov],ylab=paste(i,"->",j))
            if(plotCI){
              dN<-t(apply(desMat,1,function(x) tryCatch(numDeriv::grad(get_gamma,beta,covs=matrix(x,nrow=1),nbStates=nbStates,i=i,j=j),error=function(e) NA)))
              se<-t(apply(dN,1,function(x) tryCatch(suppressWarnings(sqrt(x%*%Sigma[gamInd,gamInd]%*%x)),error=function(e) NA)))
              if(!all(is.na(se))) {
                lci<-1/(1+exp(-(log(trMat[i,j,]/(1-trMat[i,j,]))-quantSup*(1/(trMat[i,j,]-trMat[i,j,]^2))*se)))#trMat[i,j,]-quantSup*se[i,j]
                uci<-1/(1+exp(-(log(trMat[i,j,]/(1-trMat[i,j,]))+quantSup*(1/(trMat[i,j,]-trMat[i,j,]^2))*se)))#trMat[i,j,]+quantSup*se[i,j]
                seInd <- which(!is.na(se))
                ciInd <- seInd[which(abs(uci[seInd]-lci[seInd])>max(abs(uci[seInd]-lci[seInd]))/1000)] #to aviod un-supressable warning in arrows()
                arrows(as.numeric(tempCovs[ciInd,cov]), lci[ciInd], as.numeric(tempCovs[ciInd,cov]), uci[ciInd], length=0.025, angle=90, code=3, col=gray(.5))
              }
            }
          }
        }
        if(ncol(rawCovs)>1) mtext(paste("Transition probabilities:",paste(names(rawCovs)[-cov],"=",tmpcovs[-cov],collapse=", ")),side=3,outer=TRUE,padj=2)
        else mtext("Transition probabilities:",side=3,outer=TRUE,padj=2)
        
      }
    }
  }

  #################################
  ## Plot maps colored by states ##
  #################################
  if(all(c("x","y") %in% names(m$data)) & plotTracks){
    
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
# used in the function plot.momentuHMM.
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
#  - sepStates, breaks, hist.ylim: see arguments of plot.momentuHMM.
#  - state: if sepStates, this function needs to know which state needs to be plotted.
#  - col: colors of the state-dependent density lines

plotHist <- function (gen,genDensities,dist,message,
                      sepStates,breaks="Sturges",state=NULL,hist.ylim=NULL,col=NULL,legText, cumul=TRUE)
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
    lty <- rep(1, nbStates)
    if (cumul) {
      legText <- c(legText, "Total")
      col <- c(col, "black")
      lty <- c(lty, 2)
    }
  }

  distname <- names(dist)
  
  if(dist %in% angledists){
    h <- hist(gen,plot=F,breaks=breaks) # to determine 'breaks'
    breaks <- seq(-pi,pi,length=length(h$breaks))
    
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
        if(maxdens>ymax){
          ymax <- ifelse(maxdens<2*max(h$density),maxdens,2*max(h$density))
        }
      }
    }

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
      if(cumul){
        total <- genDensities[[1]]
        for (s in 2:nbStates) total[, 2] <- total[, 2] + genDensities[[s]][, 2]
        lines(total, lwd = 2, lty = 2)
      }
      legend("topright",legText,lwd=rep(2,nbStates),col=col,bty="n")
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
        if(maxdens>ymax){
          ymax <- ifelse(maxdens<2*max(h$density),maxdens,2*max(h$density))
        }
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
      if(cumul){
        total <- genDensities[[1]]
        for (s in 2:nbStates) total[, 2] <- total[, 2] + genDensities[[s]][, 2]
        lines(total, lwd = 2, lty = 2)
      }
      legend("topright",legText,lwd=rep(2,nbStates),col=col,bty="n")
    }
  }

}

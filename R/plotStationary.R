
#' Plot stationary state probabilities
#'
#' @param model \code{\link{momentuHMM}}, \code{\link{momentuHierHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object
#' @param covs Optional data frame consisting of a single row indicating the covariate values to be used in plots.
#' If none are specified, the means of any covariates appearing in the model are used (unless covariate is a factor, in which case the first factor in the data is used).
#' @param col Vector or colors for the states (one color per state).
#' @param plotCI Logical indicating whether to include confidence intervals in plots (default: FALSE)
#' @param alpha Significance level of the confidence intervals (if \code{plotCI=TRUE}). Default: 0.95 (i.e. 95\% CIs).
#' @param ... Additional arguments passed to \code{graphics::plot}. These can currently include \code{cex.axis}, \code{cex.lab}, \code{cex.legend}, \code{cex.main}, \code{legend.pos}, and \code{lwd}. See \code{\link[graphics]{par}}. \code{legend.pos} can be a single keyword from the list ``bottomright'', ``bottom'', ``bottomleft'', ``left'', ``topleft'', ``top'', ``topright'', ``right'', and ``center''.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plotStationary(m)
#'
#' @export
#'
plotStationary <- function (model, covs = NULL, col=NULL, plotCI=FALSE, alpha=0.95, ...) {
  UseMethod("plotStationary")
}

#' @method plotStationary momentuHMM
#' @export
plotStationary.momentuHMM <- function(model, covs = NULL, col=NULL, plotCI=FALSE, alpha=0.95, ...)
{
    model <- delta_bc(model)

    data <- model$data
    nbStates <- length(model$stateNames)
    beta <- model$mle$beta

    if(nrow(beta)/model$conditions$mixtures==1)
        stop("No covariate effect to plot")
    else if(inherits(model,"hierarchical")){
      if(nrow(beta)/model$conditions$mixtures==nlevels(model$data$level)) stop("No covariate effect to plot")
    }
    
    # prepare colors for the states
    if(!is.null(col) & length(col)!=nbStates) {
        warning("Length of 'col' should be equal to number of states - argument ignored")
        col <- NULL
    }
    if(is.null(col) & nbStates<8) {
        # color-blind friendly palette
        pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        col <- pal[1:nbStates]
    }
    if(is.null(col) & nbStates>=8) {
        # to make sure that all colours are distinct (emulate ggplot default palette)
        hues <- seq(15, 375, length = nbStates + 1)
        col <- hcl(h = hues, l = 65, c = 100)[1:nbStates]
    }

    if(inherits(model,"miSum") & plotCI){
      if(length(model$conditions$optInd)){
        Sigma <- matrix(0,length(model$mod$estimate),length(model$mod$estimate))
        Sigma[(1:length(model$mod$estimate))[-model$conditions$optInd],(1:length(model$mod$estimate))[-model$conditions$optInd]] <- model$MIcombine$variance
      } else {
        Sigma <- model$MIcombine$variance
      }
    } else if(!is.null(model$mod$hessian) & plotCI){
        Sigma <- model$mod$Sigma
    } else {
        Sigma <- NULL
        plotCI <- FALSE
    }

    formula<-model$conditions$formula
    newForm <- newFormulas(formula,nbStates,model$conditions$betaRef,hierarchical=TRUE)
    newformula <- newForm$newformula
    recharge <- hierRecharge <- newForm$recharge
    
    covs <- getCovs(model,covs,unique(model$data$ID))
    
    aInd <- NULL
    nbAnimals <- length(unique(data$ID))
    for(i in 1:nbAnimals){
      aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
    }
    
    if(!is.null(recharge)){
      reForm <- formatRecharge(nbStates,formula,model$conditions$betaRef,data=data,covs=covs,par=list(g0=model$mle$g0,theta=model$mle$theta))
      newformula <- reForm$newformula
      recharge <- reForm$recharge
      hierRecharge <- reForm$hierRecharge
      data[colnames(reForm$newdata)] <- reForm$newdata
      g0covs <- reForm$g0covs
      nbG0covs <- ncol(g0covs)-1
      recovs <- reForm$recovs
      nbRecovs <- ncol(recovs)-1
      if(inherits(data,"hierarchical")) rechargeNames <- paste0("recharge",gsub("level","",names(hierRecharge)))
      else rechargeNames <- "recharge"
      
      for(j in rechargeNames){
        if(!is.null(covs[[j]])) reForm$covs[[j]] <- covs[[j]]
      }
      covs <- reForm$covs

      covsCol <- cbind(get_all_vars(newformula,data),get_all_vars(recharge$theta,data))#rownames(attr(stats::terms(formula),"factors"))#attr(stats::terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
      if(!all(names(covsCol) %in% names(data))){
        covsCol <- covsCol[,names(covsCol) %in% names(data),drop=FALSE]
      }
      rawCovs <- covsCol[,c(unique(colnames(covsCol))),drop=FALSE]
    } else {
      nbG0covs <- 0
      nbRecovs <- 0
      g0covs <- NULL
      recovs <- NULL
      rawCovs <- model$rawCovs
    }
    
    if(inherits(model,"hierarchical")) {
      covIndex <- which(!(names(rawCovs)=="level"))
      #covs$level <- NULL
      #covs <- data.frame(covs[rep(1:nrow(covs),nlevels(model$data$level)),,drop=FALSE],level=rep(levels(model$data$level),each=nrow(covs)))
    } else covIndex <- 1:ncol(rawCovs)
    
    nbCovs <- ncol(stats::model.matrix(newformula,data))-1 # substract intercept column
    mixtures <- model$conditions$mixtures
    
    gamInd<-(length(model$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)*mixtures+1):(length(model$mod$estimate))-(ncol(model$covsPi)*(mixtures-1))-ifelse(nbRecovs,(nbRecovs+1)+(nbG0covs+1),0)-ncol(model$covsDelta)*(nbStates-1)*(!model$conditions$stationary)*mixtures
    #gamInd <- gamInd[model$conditions$betaCons]
    
    # loop over covariates
    for(cov in covIndex) {

        if(!is.factor(rawCovs[,cov])){

          gridLength <- 101
          hGridLength <- gridLength#gridLength*ifelse(inherits(model,"hierarchical"),nlevels(model$data$level),1)

          inf <- min(rawCovs[,cov],na.rm=TRUE)
          sup <- max(rawCovs[,cov],na.rm=TRUE)

          # set all covariates to their mean, except for "cov"
          # (which takes a grid of values from inf to sup)
          tempCovs <- data.frame(matrix(covs[names(rawCovs)][[1]],nrow=hGridLength,ncol=1))
          if(ncol(rawCovs)>1)
            for(i in 2:ncol(rawCovs))
              tempCovs <- cbind(tempCovs,rep(covs[names(rawCovs)][[i]],gridLength))

          tempCovs[,cov] <- rep(seq(inf,sup,length=gridLength),each=hGridLength/gridLength)
        } else {
          gridLength<- nlevels(rawCovs[,cov])
          hGridLength <- gridLength#gridLength*ifelse(inherits(model,"hierarchical"),nlevels(model$data$level),1)
          # set all covariates to their mean, except for "cov"
          tempCovs <- data.frame(matrix(covs[names(rawCovs)][[1]],nrow=hGridLength,ncol=1))
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

        tmpSplineInputs<-getSplineFormula(newformula,data,tempCovs)
        
        if(inherits(model,"hierarchical")) {
          #tmpSplineInputs$covs$level <- NULL
          tempCovs <- tempCovs[rep(1:nrow(tempCovs),each=nlevels(model$data$level)),,drop=FALSE]
          tempCovs$level <- rep(levels(model$data$level),times=nrow(tmpSplineInputs$covs))
          tmpcovs <- tmpcovs[rep(1,nlevels(model$data$level)),,drop=FALSE]
          tmpcovs$level <- levels(model$data$level)
          if(is.null(recharge)){
            tmpSplineInputs$covs <- tempCovs
          }
          class(tempCovs) <- append("hierarchical",class(tempCovs))
          class(tmpSplineInputs$covs) <- append("hierarchical",class(tmpSplineInputs$covs))
        }
        statPlot(model,Sigma,nbStates,tmpSplineInputs$formula,tmpSplineInputs$covs,tempCovs,tmpcovs,cov,hierRecharge,alpha,gridLength,gamInd,names(rawCovs),col,plotCI,...)
    }
}

# for differentiation in delta method
get_stat <- function(beta,covs,nbStates,i,betaRef,betaCons,workBounds=matrix(c(-Inf,Inf),length(betaCons),2,byrow=TRUE),mixture=1,ref=1:nbStates) {
  gamma <- get_gamma(beta,covs,nbStates,1:nbStates,1:nbStates,betaRef,betaCons,workBounds,mixture)
  solve(t(diag(length(ref))-gamma[ref,ref]+1),rep(1,length(ref)))[i]
}

get_stat_recharge <- function(beta,covs,formula,hierRecharge,nbStates,i,betaRef,betaCons,workBounds=matrix(c(-Inf,Inf),length(betaCons)+length(beta[(max(betaCons)+1):length(beta)]),2,byrow=TRUE),mixture=1,ref=1:nbStates){
  gamma <- get_gamma_recharge(beta,covs,formula,hierRecharge,nbStates,1:nbStates,1:nbStates,betaRef,betaCons,workBounds,mixture)
  solve(t(diag(length(ref))-gamma[ref,ref]+1),rep(1,length(ref)))[i]
}

statPlot<-function(model,Sigma,nbStates,formula,covs,tempCovs,tmpcovs,cov,hierRecharge,alpha,gridLength,gamInd,covnames,col,plotCI,...){

    if(!missing(...)){
        arg <- list(...)
        if(any(!(names(arg) %in% plotArgs[-match(c("cex","asp"),plotArgs,nomatch=0)]))) stop("additional graphical parameters for plotStationary are currently limited to: ",paste0(plotArgs[-match(c("cex","asp"),plotArgs,nomatch=0)],collapse=", "))
        if(!is.null(arg$cex.main)) cex.main <- arg$cex.main
        else cex.main <- 1
        arg$cex.main <- NULL
        if(!is.null(arg$cex.legend)) cex.legend <- arg$cex.legend
        else cex.legend <- 1
        arg$cex.legend <- NULL
        cex <- 0.6
        arg$cex <- NULL
        asp <- 1
        arg$asp <- NULL
        if(!is.null(arg$lwd)) lwd <- arg$lwd
        else lwd <- 1.3
        arg$lwd <- NULL
        if(!is.null(arg$legend.pos)) {
          if(!(arg$legend.pos %in% c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")))
            stop('legend.pos must be a single keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"')
          legend.pos <- arg$legend.pos
        }
        else legend.pos <- NULL
        arg$legend.pos <- NULL
    } else {
        cex <- 0.6
        asp <- 1
        lwd <- 1.3
        cex.main <- 1
        cex.legend <- 1
        legend.pos <- NULL
        arg <- NULL
    }
    marg <- arg
    marg$cex <- NULL
    
    if(is.null(hierRecharge)){
      desMat <- stats::model.matrix(formula,data=covs)
      probs <- stationary(model, covs=desMat)
    } else {
      if(inherits(model,"hierarchical")) covs$level <- NULL
      probs <- stationary(model, covs=covs)
    }
    
    mixtures <- model$conditions$mixtures
    
    for(mix in 1:mixtures){
      if(!inherits(model,"hierarchical")){
        plotCall(cov,tempCovs,probs[[mix]],model,nbStates,covnames,lwd,arg,col,legend.pos,cex.legend,plotCI,gridLength,hierRecharge,desMat,mix,Sigma,gamInd,alpha,1:nbStates,model$stateNames,formula)
        if(length(covnames)>1) do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities: ",paste(covnames[-cov]," = ",tmpcovs[-cov],collapse=", ")),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
        else do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities"),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
      } else {
        for(j in 1:(model$conditions$hierStates$height-1)){
          if(j==1) {
            # only plot if there is variation in stationary state proabilities
            if(!all(apply(probs[[mix]][["level1"]],2,function(x) all( abs(x - mean(x)) < 1.e-6 )))){
              ref <- model$conditions$hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)
              plotCall(cov,tempCovs[which(tempCovs$level==j),],probs[[mix]][["level1"]],model,nbStates,covnames,lwd,arg,col,legend.pos,cex.legend,plotCI,gridLength,hierRecharge,desMat[which(tempCovs$level==j),],mix,Sigma,gamInd,alpha,ref,names(ref),formula)
              if(length(covnames[-cov][which(covnames[-cov]!="level" & !grepl("recharge",covnames[-cov]))])) do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities for level",j,": ",paste(covnames[-cov][which(covnames[-cov]!="level" & !grepl("recharge",covnames[-cov]))]," = ",tmpcovs[which(tmpcovs$level==j),-cov][which(covnames[-cov]!="level" & !grepl("recharge",covnames[-cov]))],collapse=", ")),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
              else do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities for level",j),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
            }
          } else {
            t <- data.tree::Traverse(model$conditions$hierStates,filterFun=function(x) x$level==j)
            names(t) <- model$conditions$hierStates$Get("name",filterFun=function(x) x$level==j)
            for(k in names(t)){
              ref <- t[[k]]$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)#t[[k]]$Get("state",filterFun = data.tree::isLeaf)
              # only plot if jth node has children and there is variation in stationary state proabilities
              if(!is.null(ref) && !all(apply(probs[[mix]][[paste0("level",j)]][[k]],2,function(x) all( abs(x - mean(x)) < 1.e-6 )))){
                plotCall(cov,tempCovs[which(tempCovs$level==j),],probs[[mix]][[paste0("level",j)]][[k]],model,nbStates,covnames,lwd,arg,col,legend.pos,cex.legend,plotCI,gridLength,hierRecharge,desMat[which(tempCovs$level==j),],mix,Sigma,gamInd,alpha,ref,names(ref),formula)
                if(length(covnames[-cov][which(covnames[-cov]!="level" & !grepl("recharge",covnames[-cov]))])) do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities for level",j," ",k,": ",paste(covnames[-cov][which(covnames[-cov]!="level" & !grepl("recharge",covnames[-cov]))]," = ",tmpcovs[which(tmpcovs$level==j),-cov][which(covnames[-cov]!="level" & !grepl("recharge",covnames[-cov]))],collapse=", ")),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
                else do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities for level",j," ",k),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
              }
            }
          }
        }
      }
    }
}

plotCall <- function(cov,tempCovs,pr,model,nbStates,covnames,lwd,arg,col,legend.pos,cex.legend,plotCI,gridLength,hierRecharge,desMat,mix,Sigma,gamInd,alpha,ref=1:nbStates,stateNames,formula){
  if(!is.factor(tempCovs[,cov])){
    do.call(plot,c(list(tempCovs[,cov],pr[,1],type="l",ylim=c(0,1),col=col[ref[1]],xlab=covnames[cov], ylab="Stationary state probabilities",lwd=lwd),arg))
    for(state in 2:length(ref))
      points(tempCovs[,cov], pr[,state], type="l", col=col[ref[state]])
  } else {
    do.call(plot,c(list(tempCovs[,cov],pr[,1],type="l",ylim=c(0,1),col=col[ref[1]],xlab=covnames[cov], ylab="Stationary state probabilities",lwd=lwd,border=col[ref[1]]),arg))
    for(state in 2:length(ref))
      plot(tempCovs[,cov], pr[,state], type="l", col=col[ref[state]], border=col[ref[state]], add=TRUE)
  }
  legend(ifelse(is.null(legend.pos),"topleft",legend.pos),stateNames,lwd=rep(lwd,length(ref)),col=col[ref],lty=1,bty="n",cex=cex.legend)
  
  if(plotCI) {
    
    # this function is used to muffle the warning "zero-length arrow is of indeterminate angle and so skipped" when plotCI=TRUE
    muffWarn <- function(w) {
      if(any(grepl("zero-length arrow is of indeterminate angle and so skipped",w)))
        invokeRestart("muffleWarning")
    }
    
    lci <- matrix(NA,gridLength,length(ref))
    uci <- matrix(NA,gridLength,length(ref))
    
    for(state in 1:length(ref)) {
      
      if(is.null(hierRecharge)){
        dN <- t(apply(desMat, 1, function(x)
          numDeriv::grad(get_stat,model$mod$estimate[gamInd[unique(c(model$conditions$betaCons))]],covs=matrix(x,1),nbStates=nbStates,i=state,betaRef=model$conditions$betaRef,betaCons=model$conditions$betaCons,workBounds=model$conditions$workBounds$beta,mixture=mix,ref=ref)))
        tmpSig <- Sigma[gamInd[unique(c(model$conditions$betaCons))],gamInd[unique(c(model$conditions$betaCons))]]
      } else {
        recharge <- expandRechargeFormulas(hierRecharge)
        recovs <- stats::model.matrix(recharge$theta,tempCovs)
        nbRecovs <- ncol(recovs)-1
        tmpSig <- Sigma[c(gamInd[unique(c(model$conditions$betaCons))],length(model$mod$estimate)-nbRecovs:0),c(gamInd[unique(c(model$conditions$betaCons))],length(model$mod$estimate)-nbRecovs:0)]
        dN<-matrix(unlist(lapply(split(tempCovs,1:nrow(tempCovs)),function(x) tryCatch(numDeriv::grad(get_stat_recharge,model$mod$estimate[c(gamInd[unique(c(model$conditions$betaCons))],length(model$mod$estimate)-nbRecovs:0)],covs=x,formula=formula,hierRecharge=hierRecharge,nbStates=nbStates,i=state,betaRef=model$conditions$betaRef,betaCons=model$conditions$betaCons,workBounds=rbind(model$conditions$workBounds$beta,model$conditions$workBounds$theta),mixture=mix,ref=ref),error=function(e) NA))),ncol=ncol(tmpSig),byrow=TRUE)
      }
      
      se <- t(apply(dN, 1, function(x)
        suppressWarnings(sqrt(x%*%tmpSig%*%x))))
      
      lci[,state] <- 1/(1 + exp(-(log(pr[,state]/(1-pr[,state])) -
                                    qnorm(1-(1-alpha)/2) * (1/(pr[,state]-pr[,state]^2)) * se)))
      uci[,state] <- 1/(1 + exp(-(log(pr[,state]/(1-pr[,state])) +
                                    qnorm(1-(1-alpha)/2) * (1/(pr[,state]-pr[,state]^2)) * se)))
      
      # plot the confidence intervals
      ciInd <- which(!is.na(se))
      
      withCallingHandlers(do.call(arrows,c(list(as.numeric(tempCovs[ciInd,cov]), lci[ciInd,state], as.numeric(tempCovs[ciInd,cov]),
                                                uci[ciInd,state], length=0.025, angle=90, code=3, col=col[ref[state]], lwd=lwd),arg)),warning=muffWarn)
      
    }
  }
}

#' @method plotStationary miSum
#' @export
plotStationary.miSum <- function(model, covs = NULL, col=NULL, plotCI=FALSE, alpha=0.95, ...)
{
  model <- delta_bc(model)
  model <- formatmiSum(model)
  plotStationary(momentuHMM(model),covs,col,plotCI,alpha,...)
}

#' @method plotStationary miHMM
#' @export
plotStationary.miHMM <- function(model, covs = NULL, col=NULL, plotCI=FALSE, alpha=0.95, ...)
{
  plotStationary(model$miSum,covs,col,plotCI,alpha,...)
}

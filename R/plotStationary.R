
#' Plot stationary state probabilities
#'
#' @param model \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object
#' @param covs Optional data frame consisting of a single row indicating the covariate values to be used in plots.
#' If none are specified, the means of any covariates appearing in the model are used (unless covariate is a factor, in which case the first factor in the data is used).
#' @param col Vector or colors for the states (one color per state).
#' @param plotCI Logical indicating whether to include confidence intervals in plots (default: FALSE)
#' @param alpha Significance level of the confidence intervals (if \code{plotCI=TRUE}). Default: 0.95 (i.e. 95\% CIs).
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}. These can currently include \code{cex.axis}, \code{cex.lab}, \code{cex.legend}, \code{cex.main}, \code{legend.pos}, and \code{lwd}. See \code{\link[graphics]{par}}. \code{legend.pos} can be a single keyword from the list ``bottomright'', ``bottom'', ``bottomleft'', ``left'', ``topleft'', ``top'', ``topright'', ``right'', and ``center''.
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

    if(nrow(beta)==1)
        stop("No covariate effect to plot (nrow(beta)==1).")

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

    rawCovs <- model$rawCovs

    if(inherits(model,"miSum") & plotCI){
        Sigma <- model$MIcombine$variance
    } else if(!is.null(model$mod$hessian) & plotCI){
        Sigma <- ginv(model$mod$hessian)
    } else {
        Sigma <- NULL
        plotCI <- FALSE
    }

    # indices corresponding to regression coefficients in model$mod$estimate
    formula<-model$conditions$formula
    newForm <- newFormulas(formula,nbStates)
    formulaStates <- newForm$formulaStates
    formterms <- newForm$formterms
    newformula <- newForm$newformula
    
    nbCovs <- ncol(model.matrix(newformula,model$data))-1 # substract intercept column
    gamInd<-(length(model$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)+1):(length(model$mod$estimate))-ncol(model$covsDelta)*(nbStates-1)*(!model$conditions$stationary)

    if(is.null(covs)){
        covs <- model$data[1,]
        for(j in names(model$data)[which(unlist(lapply(model$data,function(x) any(class(x) %in% meansList))))]){
          if(inherits(model$data[[j]],"angle")) covs[[j]] <- CircStats::circ.mean(model$data[[j]][!is.na(model$data[[j]])])
          else covs[[j]]<-mean(model$data[[j]],na.rm=TRUE)
        }
    } else {
        if(!is.data.frame(covs)) stop('covs must be a data frame')
        if(nrow(covs)>1) stop('covs must consist of a single row')
        if(!all(names(covs) %in% names(model$data))) stop('invalid covs specified')
        if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(model$data$ID))
        for(j in names(model$data)[which(names(model$data) %in% names(covs))]){
          if(inherits(model$data[[j]],"factor")) covs[[j]] <- factor(covs[[j]],levels=levels(model$data[[j]]))
          if(is.na(covs[[j]])) stop("check value for ",j)
        }
        for(j in names(model$data)[which(!(names(model$data) %in% names(covs)))]){
          if(inherits(model$data[[j]],"factor")) covs[[j]] <- model$data[[j]][1]
          else if(inherits(model$data[[j]],"angle")) covs[[j]] <- CircStats::circ.mean(model$data[[j]][!is.na(model$data[[j]])])
          else if(any(class(model$data[[j]]) %in% meansList)) covs[[j]]<-mean(model$data[[j]],na.rm=TRUE)
        }
    }

    # loop over covariates
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

        tmpSplineInputs<-getSplineFormula(newformula,model$data,tempCovs)

        desMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)

        statPlot(model,beta,Sigma,nbStates,desMat,tempCovs,tmpcovs,cov,alpha,gridLength,gamInd,names(rawCovs),col,plotCI,...)
    }
}

# for differentiation in delta method
get_stat <- function(beta,covs,nbStates,i) {
    gamma <- trMatrix_rcpp(nbStates,beta,covs)[,,1]
    solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))[i]
}

statPlot<-function(model,beta,Sigma,nbStates,desMat,tempCovs,tmpcovs,cov,alpha,gridLength,gamInd,covnames,col,plotCI,...){

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

    probs <- stationary(model, covs=desMat)

    if(!is.factor(tempCovs[,cov])){
        do.call(plot,c(list(tempCovs[,cov],probs[,1],type="l",ylim=c(0,1),col=col[1],xlab=covnames[cov], ylab="Stationary state probabilities",lwd=lwd),arg))
        for(state in 2:nbStates)
            points(tempCovs[,cov], probs[,state], type="l", col=col[state])
    } else {
        do.call(plot,c(list(tempCovs[,cov],probs[,1],type="l",ylim=c(0,1),col=col[1],xlab=covnames[cov], ylab="Stationary state probabilities",lwd=lwd,border=col[1]),arg))
        for(state in 2:nbStates)
            plot(tempCovs[,cov], probs[,state], type="l", col=col[state], border=col[state], add=TRUE)
    }
    legend(ifelse(is.null(legend.pos),"topleft",legend.pos),model$stateNames,lwd=rep(lwd,nbStates),col=col,lty=1,bty="n",cex=cex.legend)

    if(plotCI) {

        # this function is used to muffle the warning "zero-length arrow is of indeterminate angle and so skipped" when plotCI=TRUE
        muffWarn <- function(w) {
           if(any(grepl("zero-length arrow is of indeterminate angle and so skipped",w)))
           invokeRestart("muffleWarning")
        }

        lci <- matrix(NA,gridLength,nbStates)
        uci <- matrix(NA,gridLength,nbStates)

        for(state in 1:nbStates) {
            dN <- t(apply(desMat, 1, function(x)
                numDeriv::grad(get_stat,beta,covs=matrix(x,nrow=1),nbStates=nbStates,i=state)))

            se <- t(apply(dN, 1, function(x)
                suppressWarnings(sqrt(x%*%Sigma[gamInd,gamInd]%*%x))))

            lci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) -
                                            qnorm(1-(1-alpha)/2) * (1/(probs[,state]-probs[,state]^2)) * se)))
            uci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) +
                                            qnorm(1-(1-alpha)/2) * (1/(probs[,state]-probs[,state]^2)) * se)))

            # plot the confidence intervals
            ciInd <- which(!is.na(se))

            withCallingHandlers(do.call(arrows,c(list(as.numeric(tempCovs[ciInd,cov]), lci[ciInd,state], as.numeric(tempCovs[ciInd,cov]),
                                           uci[ciInd,state], length=0.025, angle=90, code=3, col=col[state], lwd=lwd),arg)),warning=muffWarn)

        }
    }
    if(length(covnames)>1) do.call(mtext,c(list(paste("Stationary state probabilities:",paste(covnames[-cov],"=",tmpcovs[-cov],collapse=", ")),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
    else do.call(mtext,c(list("Stationary state probabilities",side=3,outer=TRUE,padj=2,cex=cex.main),marg))
}

#' @method plotStationary miSum
#' @export
plotStationary.miSum <- function(model, covs = NULL, col=NULL, plotCI=FALSE, alpha=0.95, ...)
{
  model <- delta_bc(model)
  model$mle <- lapply(model$Par$real,function(x) x$est)
  model$mle$beta <- model$Par$beta$beta$est
  model$mle$delta <- model$Par$real$delta$est
  model$mod <- list()
  model$mod$estimate <- model$MIcombine$coefficients
  plotStationary(momentuHMM(model),covs,col,plotCI,alpha,...)
}

#' @method plotStationary miHMM
#' @export
plotStationary.miHMM <- function(model, covs = NULL, col=NULL, plotCI=FALSE, alpha=0.95, ...)
{
  plotStationary(model$miSum,covs,col,plotCI,alpha,...)
}

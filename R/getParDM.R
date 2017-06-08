#' Get starting values on working scale based on design matrix and other parameter constraints
#' 
#' Convert starting values on the natural scale of data stream probability distributions to
#' a feasible set of working scale parameters based on a design matrix and other parameter constraints.
#' 
#' If design matrix includes non-factor covariates, then natural scale parameters are assumed to correspond to the 
#' mean value(s) for the covariate(s) (if \code{nrow(data)>1}) and \code{getParDM} simply returns one possible solution to the 
#' system of linear equations defined by \code{Par}, \code{DM}, and any other constraints using singular value decomposition. 
#' This can be helpful for exploring relationships between the natural and working scale parameters when covariates are included, but \code{getParDM}
#' will not necessarily return ``good'' starting values (i.e., \code{Par0}) for \code{\link{fitHMM}} or \code{\link{MIfitHMM}}. 
#'
#' @param data Optional \code{\link{momentuHMMData}} object or a data frame containing the covariate values. 
#' \code{data} must be specified if covariates are included in \code{DM}.
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. Currently
#' supported distributions are 'gamma','weibull','exp','lnorm','beta','pois','wrpcauchy', and 'vm'. For example,
#' \code{dist=list(step='gamma', angle='vm', dives='pois')} indicates 3 data streams ('step', 'angle', and 'dives')
#' and their respective probability distributions ('gamma', 'vm', and 'pois').
#' @param Par A named list containing vectors of state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. The parameters should be on the natural scale,
#' in the order expected by the pdfs of \code{dist}, and any zero-mass parameters should be the last.
#' @param zeroInflation A named list of logicals indicating whether the probability distributions of the data streams should be zero-inflated. If \code{zeroInflation} is \code{TRUE} 
#' for a given data stream, then values for the zero-mass parameters should be
#' included in the corresponding element of \code{Par}. Ignored if \code{data} is a \code{\link{momentuHMMData}} object.
#' @param oneInflation Named list of logicals indicating whether the probability distributions of the data streams are one-inflated. If \code{oneInflation} is \code{TRUE} 
#' for a given data stream, then values for the one-mass parameters should be
#' included in the corresponding element of \code{Par}. Ignored if \code{data} is a \code{\link{momentuHMMData}} object.
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). Any \code{estAngleMean} elements corresponding to data streams that do not have angular distributions are ignored.
#' @param circularAngleMean An optional named list indicating whether to use circular-linear (FALSE) or circular-circular (TRUE) 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles. \code{circularAngleMean} elements corresponding to angular data 
#' streams are ignored unless the corresponding element of \code{estAngleMean} is \code{TRUE}. Any \code{circularAngleMean} elements 
#' corresponding to data streams that do not have angular distributions are ignored.
#' @param DM A named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of linear regression formulas or a matrix.  For example, for a 2-state 
#' model using the gamma distribution for a data stream named 'step', \code{DM=list(step=list(mean=~cov1, sd=~1))} specifies the mean 
#' parameters as a function of the covariate 'cov1' for each state.  This model could equivalently be specified as a 4x6 matrix using 
#' character strings for the covariate: 
#' \code{DM=list(step=matrix(c(1,0,0,0,'cov1',0,0,0,0,1,0,0,0,'cov1',0,0,0,0,1,0,0,0,0,1),4,6))}
#' where the 4 rows correspond to the state-dependent paramaters (mean_1,mean_2,sd_1,sd_2) and the 6 columns correspond to the regression 
#' coefficients. 
#' @param cons An optional named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. While there could be other uses, primarily intended to constrain specific parameters to be positive. For example, 
#' \code{cons=list(step=c(1,2,1,1))} raises the second parameter to the second power. Default=NULL, which simply raises all parameters to 
#' the power of 1. \code{cons} is ignored for any given data stream unless \code{DM} is specified.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. For example, for a 2-state model using the wrapped Cauchy ('wrpcauchy') distribution for 
#' a data stream named 'angle' with \code{estAngleMean$angle=TRUE)}, \code{userBounds=list(angle=matrix(c(-pi,-pi,-1,-1,pi,pi,1,1),4,2))} 
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param workcons An optional named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. Warning: use of \code{workcons} is recommended only for advanced users implementing unusual parameter constraints 
#' through a combination of \code{DM}, \code{cons}, and \code{workcons}. \code{workcons} is ignored for any given data stream unless \code{DM} is specified.
#'
#' @return A list of parameter values that can be used as starting values (\code{Par0}) in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}
#' 
#' @seealso \code{\link{getPar0}}, \code{\link{fitHMM}}, \code{\link{MIfitHMM}}
#'
#' @examples
#' # data is a momentuHMMData object, automatically loaded with the package
#' data <- example$m$data
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' nbStates <- 2
#' stepPar0 <- c(15,50,10,20) # natural scale mean_1, mean_2, sd_1, sd_2
#' anglePar0 <- c(0.7,1.5) # natural scale conentration_1, concentration_2
#' 
#' # get working parameters for 'DM' and 'cons' that constrain step length mean_1 < mean_2
#' stepDM <- matrix(c(1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,
#'           dimnames=list(NULL,c("mean:(Intercept)","mean_2",
#'                                "sd_1:(Intercept)","sd_2:(Intercept)")))
#' stepcons <- c(1,2,1,1) # coefficient for 'mean_2' constrained to be positive
#' wPar0 <- getParDM(nbStates=2,dist=list(step=stepDist),
#'                       Par=list(step=stepPar0),
#'                       DM=list(step=stepDM),cons=list(step=stepcons))
#'
#' \dontrun{
#' # Fit HMM using wPar0 as initial values for the step data stream
#' mPar <- fitHMM(data,nbStates=2,dist=list(step=stepDist,angle=angleDist),
#'                Par0=list(step=wPar0$step,angle=anglePar0),
#'                DM=list(step=stepDM),cons=list(step=stepcons))
#' }
#' 
#' # get working parameters for 'DM' using 'cov1' and 'cov2' covariates
#' stepDM2 <- list(mean=~cov1,sd=~cov2)
#' wPar20 <- getParDM(data,nbStates=2,dist=list(step=stepDist),
#'                       Par=list(step=stepPar0),
#'                       DM=list(step=stepDM2))
#'
#' \dontrun{
#' # Fit HMM using wPar20 as initial values for the step data stream
#' mPar2 <- fitHMM(data,nbStates=2,dist=list(step=stepDist,angle=angleDist),
#'                Par0=list(step=wPar20$step,angle=anglePar0),
#'                DM=list(step=stepDM2))
#' }
#'
#' @importFrom CircStats circ.mean
#' @importFrom nleqslv nleqslv
#' @export
getParDM<-function(data=data.frame(),nbStates,dist,
                 Par,
                 zeroInflation=NULL,
                 oneInflation=NULL,
                 estAngleMean=NULL,
                 circularAngleMean=NULL,
                 DM=NULL,cons=NULL,userBounds=NULL,workcons=NULL){
  
  ## check that the data is a momentuHMMData object or valid data frame
  if(!is.momentuHMMData(data))
    if(!is.data.frame(data)) stop('data must be a data.frame')
    #else if(ncol(data)<1) stop('data is empty')
  
  if(nbStates<1) stop('nbStates must be >0')
  
  if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
  if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
  distnames<-names(dist)
  if(!all(distnames %in% names(Par))) stop(distnames[which(!(distnames %in% names(Par)))]," is missing in 'Par'")
  Par <- Par[distnames]
  
  if(is.momentuHMMData(data)){
    if(any(is.na(match(distnames,names(data))))) stop(paste0(distnames[is.na(match(distnames,names(data)))],collapse=", ")," not found in data")
    
    zeroInflation <- vector('list',length(distnames))
    names(zeroInflation) <- distnames
    for(i in distnames){
      if(dist[[i]] %in% zeroInflationdists){
        if(length(which(data[[i]]==0))>0) {
          zeroInflation[[i]]<-TRUE
        }
        else 
          zeroInflation[[i]]<-FALSE
      }
      else zeroInflation[[i]]<-FALSE
    }
    oneInflation <- vector('list',length(distnames))
    names(oneInflation) <- distnames
    for(i in distnames){
      if(dist[[i]] %in% oneInflationdists){
        if(length(which(data[[i]]==1))>0) {
          oneInflation[[i]]<-TRUE
        }
        else 
          oneInflation[[i]]<-FALSE
      }
      else oneInflation[[i]]<-FALSE
    }
  } else {
    if(is.null(zeroInflation)){
      zeroInflation <- vector('list',length(distnames))
      names(zeroInflation) <- distnames
      for(i in distnames){
        zeroInflation[[i]]<-FALSE
      }
    } else {
      if(!is.list(zeroInflation) | is.null(names(zeroInflation))) stop("'zeroInflation' must be a named list")
      for(i in distnames){
        if(is.null(zeroInflation[[i]])) zeroInflation[[i]] <- FALSE
      }
    }
    if(is.null(oneInflation)){
      oneInflation <- vector('list',length(distnames))
      names(oneInflation) <- distnames
      for(i in distnames){
        oneInflation[[i]]<-FALSE
      }
    } else {
      if(!is.list(oneInflation) | is.null(names(oneInflation))) stop("'oneInflation' must be a named list")
      for(i in distnames){
        if(is.null(oneInflation[[i]])) oneInflation[[i]] <- FALSE
      }
    }
    
    if(!all(unlist(lapply(zeroInflation,is.logical)))) stop("zeroInflation must be a list of logical objects")
    if(!all(unlist(lapply(oneInflation,is.logical)))) stop("oneInflation must be a list of logical objects")
    for(i in distnames){
      if(!(dist[[i]] %in% zeroInflationdists) & zeroInflation[[i]])
        stop(dist[[i]]," distribution cannot be zero inflated")
      if(!(dist[[i]] %in% oneInflationdists) & oneInflation[[i]])
        stop(dist[[i]]," distribution cannot be one inflated")
    }
  }
  
  tempCovs <- data[1,,drop=FALSE]
  if(length(data)){
    for(j in names(data)[which(unlist(lapply(data,function(x) any(class(x) %in% meansList))))]){
      if(inherits(data[[j]],"angle")) tempCovs[[j]] <- CircStats::circ.mean(data[[j]][!is.na(data[[j]])])
      else tempCovs[[j]]<-mean(data[[j]],na.rm=TRUE)
    }
  }
  inputs <- checkInputs(nbStates,dist,Par,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,cons,workcons,stateNames=NULL)
  
  DMinputs<-getDM(tempCovs,inputs$DM,dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,inputs$cons,inputs$workcons,zeroInflation,oneInflation,inputs$circularAngleMean,FALSE)
  fullDM <- DMinputs$fullDM
  if(length(data))
    DMind <- getDM(data,inputs$DM,dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,inputs$cons,inputs$workcons,zeroInflation,oneInflation,inputs$circularAngleMean,FALSE)$DMind
  else DMind <- DMinputs$DMind
    cons <- DMinputs$cons
  workcons <- DMinputs$workcons
  
  wpar <- Par
  for(i in distnames){
    if(!is.null(DM[[i]])){
      if(DMind[[i]] & nrow(unique(fullDM[[i]]))==ncol(unique(fullDM[[i]])) & !inputs$circularAngleMean[[i]]){
        if(length(wpar[[i]])!=nrow(fullDM[[i]])) stop('Par$',i,' should be of length ',nrow(fullDM[[i]]))
        bounds<-inputs$p$bounds[[i]]
        if(any(wpar[[i]]<=bounds[,1] | wpar[[i]]>=bounds[,2])) stop('Par$',i,' must be within parameter bounds')
        bndInd <- which(!duplicated(getboundInd(fullDM[[i]])))
        a<-bounds[bndInd,1]
        b<-bounds[bndInd,2]
        par <- wpar[[i]][bndInd]
        if(any(wpar[[i]]!=par[getboundInd(fullDM[[i]])])) stop('Par$',i,' values are not consistent with DM$',i)
        piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
        ind1<-which(piInd)
        ind2<-which(!piInd)
        
        p<-numeric(length(bndInd))
        if(length(ind1)){
          if(!inputs$circularAngleMean[[i]]) p[ind1] <- (solve(unique(fullDM[[i]])[ind1,ind1],tan(par[ind1]/2))-workcons[[i]][ind1])^(1/cons[[i]][ind1])
          else stop("sorry, circular angle mean parameters are not supported by getParDM")
        }
        
        ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
        ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
        ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
        
        if(length(ind21)) p[ind21]<-(solve(unique(fullDM[[i]])[ind21,ind21],log(par[ind21]-a[ind21]))-workcons[[i]][ind21])^(1/cons[[i]][ind21])
        if(length(ind22)) p[ind22]<-(solve(unique(fullDM[[i]])[ind22,ind22],logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22])))-workcons[[i]][ind22])^(1/cons[[i]][ind22])
        if(length(ind23)) p[ind23]<-(solve(unique(fullDM[[i]])[ind23,ind23],-log(-par[ind23]+b[ind23]))-workcons[[i]][ind23])^(1/cons[[i]][ind23])
        
      } else {
        
        if(length(wpar[[i]])!=nrow(fullDM[[i]])) stop('Par$',i,' should be of length ',nrow(fullDM[[i]]))
        bounds<-inputs$p$bounds[[i]]
        if(any(wpar[[i]]<=bounds[,1] | wpar[[i]]>=bounds[,2])) stop('Par$',i,' must be within parameter bounds')
        if(is.list(inputs$DM[[i]])){
          if(!inputs$circularAngleMean[[i]])
            gbInd <- getboundInd(fullDM[[i]])
          else
            gbInd <- c(1:nbStates,getboundInd(fullDM[[i]][1:nbStates+nbStates,])+nbStates)
        } else {
          if(!inputs$circularAngleMean[[i]])
            gbInd <- getboundInd(inputs$DM[[i]])
          else 
            gbInd <- c(1:nbStates,getboundInd(inputs$DM[[i]][1:nbStates+nbStates,])+nbStates)
        }
        bndInd <- which(!duplicated(gbInd))
        a<-bounds[bndInd,1]
        b<-bounds[bndInd,2]
        par <- wpar[[i]][bndInd]
        if(any(wpar[[i]]!=par[gbInd])) stop('Par$',i,' values are not consistent with DM$',i)
        piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
        ind1<-which(piInd)
        ind2<-which(!piInd)
        
        ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
        ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
        ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
        
        if(length(ind1)){
          if(inputs$estAngleMean[[i]]){
            
            p<-numeric(ncol(fullDM[[i]]))
            meanind<-which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],2,function(x) !all(unlist(x)==0))))
            
            if(!inputs$circularAngleMean[[i]]) {
              
              asvd<-svd(fullDM[[i]][gbInd,][ind1,meanind])
              adiag <- diag(1/asvd$d)
              p[meanind] <- ((asvd$v %*% adiag %*% t(asvd$u) %*% tan(par[ind1]/2))-workcons[[i]][meanind])^(1/cons[[i]][meanind])
              #p[meanind] <- (solve(unique(fullDM[[i]])[ind1,meanind],tan(par[ind1]/2))-workcons[[i]][meanind])^(1/cons[[i]][meanind])

              if(length(ind21)){
                asvd<-svd(fullDM[[i]][gbInd,][ind21,-meanind])
                adiag <- diag(1/asvd$d)
                p[-meanind] <- ((asvd$v %*% adiag %*% t(asvd$u) %*% log(par[ind21]-a[ind21]))-workcons[[i]][-meanind])^(1/cons[[i]][-meanind])
              }
              
              if(length(ind22)){
                asvd<-svd(fullDM[[i]][gbInd,][ind22,-meanind])
                adiag <- diag(1/asvd$d)
                p[-meanind] <- ((asvd$v %*% adiag %*% t(asvd$u) %*% logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22])))-workcons[[i]][-meanind])^(1/cons[[i]][-meanind])
              }
              
              if(length(ind23)){
                asvd<-svd(fullDM[[i]][gbInd,][ind23,-meanind])
                adiag <- diag(1/asvd$d)
                p[-meanind] <- ((asvd$v %*% adiag %*% t(asvd$u) %*% (-log(-par[ind23]+b[ind23])))-workcons[[i]][-meanind])^(1/cons[[i]][-meanind])
              }
            } else {
              
              meanind1<-which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
              meanind2<-which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],2,function(x) !all(unlist(x)==0))))
              xmat <- fullDM[[i]][gbInd,][meanind1,meanind2]
              nc<-apply(xmat,1:2,function(x) !all(unlist(x)==0))
              
              solveatan2<-function(x,theta,covs,cons,workcons){
                Xvec <- x^cons+workcons
                XB<-rep(0,length(meanind1))
                XB1<-XB2<-rep(0,length(meanind1))
                for(i in 1:length(meanind1)){
                  ncind <- which(nc[meanind1[i],])
                  XB1[i]<-sin(covs[i,ncind])%*%Xvec[ncind]
                  XB2[i]<-cos(covs[i,ncind])%*%Xvec[ncind]
                  XB[i] <- atan2(XB1[i],1+XB2[i])
                }
                c(abs(theta - XB),rep(0,length(x)-length(theta)))
              }

              if(length(meanind1)) p[meanind2] <- nleqslv::nleqslv(x=rep(1,ncol(xmat)),fn=solveatan2,theta=par[meanind1],covs=xmat,cons=cons[[i]][meanind2],workcons=workcons[[i]][meanind2],control=list(allowSingular=TRUE))$x
              
              meanind<-which((apply(fullDM[[i]][nbStates+1:nbStates,,drop=FALSE],2,function(x) !all(unlist(x)==0))))
              
              if(length(ind21)){
                asvd<-svd(fullDM[[i]][gbInd,][ind21,meanind])
                adiag <- diag(1/asvd$d)
                p[meanind] <- ((asvd$v %*% adiag %*% t(asvd$u) %*% log(par[ind21]-a[ind21]))-workcons[[i]][meanind])^(1/cons[[i]][meanind])
              }
              
              if(length(ind22)){
                asvd<-svd(fullDM[[i]][gbInd,][ind22,meanind])
                adiag <- diag(1/asvd$d)
                p[meanind] <- ((asvd$v %*% adiag %*% t(asvd$u) %*% logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22])))-workcons[[i]][meanind])^(1/cons[[i]][meanind])
              }
              
              if(length(ind23)){
                asvd<-svd(fullDM[[i]][gbInd,][ind23,meanind])
                adiag <- diag(1/asvd$d)
                p[meanind] <- ((asvd$v %*% adiag %*% t(asvd$u) %*% (-log(-par[ind23]+b[ind23])))-workcons[[i]][meanind])^(1/cons[[i]][meanind])
              }
            
            }
          } else if(!length(ind2)){
            p <- (solve(fullDM[[i]][gbInd,],tan(par/2))-workcons[[i]])^(1/cons[[i]])
          } else stop("sorry, the parameters for ",i," cannot have different bounds")
        } else if(((length(ind21)>0) + (length(ind22)>0) + (length(ind23)>0))>1){ 
          stop("sorry, getParDM requires the parameters for ",i," to have identical bounds when covariates are included in the design matrix")
        } else {
          if(length(ind21)){
            asvd<-svd(fullDM[[i]][gbInd,][ind21,])
            adiag <- diag(1/asvd$d)
            p <- ((asvd$v %*% adiag %*% t(asvd$u) %*% log(par[ind21]-a[ind21]))-workcons[[i]])^(1/cons[[i]])
          }
          
          if(length(ind22)){
            asvd<-svd(fullDM[[i]][gbInd,][ind22,])
            adiag <- diag(1/asvd$d)
            p <- ((asvd$v %*% adiag %*% t(asvd$u) %*% logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22])))-workcons[[i]])^(1/cons[[i]])
          }
          
          if(length(ind23)){
            asvd<-svd(fullDM[[i]][gbInd,][ind23,])
            adiag <- diag(1/asvd$d)
            p <- ((asvd$v %*% adiag %*% t(asvd$u) %*% (-log(-par[ind23]+b[ind23])))-workcons[[i]])^(1/cons[[i]])
          }
        }
      }
      if(any(!is.finite(p))) stop(i," working scale parameters are not finite. Check natural parameter values, bounds, and constraints.")
      wpar[[i]]<-c(p)      
    }
  }
  wpar
}


#' Confidence intervals for the natural (i.e., real) parameters
#'
#' Computes the standard errors and confidence intervals on the real (i.e., natural) scale of the data stream probability distribution parameters,
#' as well as for the transition probabilities parameters. If covariates are included in the probability distributions or TPM formula, the mean values
#' of non-factor covariates are used for calculating the natural parameters. For any covariate(s) of class 'factor', then the value(s) from the first observation 
#' in the data are used.
#'
#' @param m A \code{momentuHMM} object
#' @param alpha Significance level of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param covs Data frame consisting of a single row indicating the covariate values to be used in the calculations. 
#' For any covariates that are not specified using \code{covs}, the means of the covariate(s) are used 
#' (unless the covariate is a factor, in which case the first factor in the data is used). By default, no covariates are specified.
#'
#' @return A list of the following objects:
#' \item{...}{List(s) of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the natural parameters of the data streams}
#' \item{gamma}{List of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the transition probabilities}
#' \item{delta}{List of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the initial state probabilities}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' ci1<-CIreal(m)
#' 
#' # specify 'covs'
#' ci2<-CIreal(m,covs=data.frame(cov1=mean(m$data$cov1),cov2=mean(m$data$cov2)))
#' 
#' all.equal(ci1,ci2)
#'
#' @export
#' @importFrom numDeriv grad
#' @importFrom utils tail
#' @importFrom Brobdingnag as.brob sum
#' @importFrom CircStats circ.mean

CIreal <- function(m,alpha=0.95,covs=NULL)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  if(length(m$mod)<=1)
    stop("The given model hasn't been fitted.")

  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1.")

  nbStates <- length(m$stateNames)

  dist <- m$conditions$dist
  distnames <- names(dist)
  DMind <- m$conditions$DMind
  
  m <- delta_bc(m)
  
  # identify covariates
  if(is.null(covs)){
    tempCovs <- m$data[1,]
    for(j in names(m$data)[which(unlist(lapply(m$data,function(x) any(class(x) %in% meansList))))]){
      if(inherits(m$data[[j]],"angle")) tempCovs[[j]] <- CircStats::circ.mean(m$data[[j]][!is.na(m$data[[j]])])
      else tempCovs[[j]]<-mean(m$data[[j]],na.rm=TRUE)
    }
    covs <- tempCovs
  } else {
    if(!is.data.frame(covs)) stop('covs must be a data frame')
    if(nrow(covs)>1) stop('covs must consist of a single row')
    if(is.null(recharge))
      if(!all(names(covs) %in% names(m$data))) stop('invalid covs specified')
    else 
      if(!all(names(covs) %in% c(names(m$data),"recharge"))) stop('invalid covs specified')
    if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(m$data$ID))
    for(j in names(m$data)[which(!(names(m$data) %in% names(covs)))]){
      if(any(class(m$data[[j]]) %in% meansList)){
        if(inherits(m$data[[j]],"angle")) covs[[j]] <- CircStats::circ.mean(m$data[[j]][!is.na(m$data[[j]])])
        else covs[[j]]<-mean(m$data[[j]],na.rm=TRUE)
      } else covs[[j]] <- m$data[[j]][1]
    }
    for(j in names(m$data)[which(names(m$data) %in% names(covs))]){
      if(inherits(m$data[[j]],"factor")) covs[[j]] <- factor(covs[[j]],levels=levels(m$data[[j]]))
      if(is.na(covs[[j]])) stop("check covs value for ",j)
    }    
    tempCovs <- covs[1,]
  }
  
  formula<-m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  nbCovs <- ncol(model.matrix(newformula,m$data))-1 # substract intercept column
  
  aInd <- NULL
  nbAnimals <- length(unique(m$data$ID))
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(m$data$ID==unique(m$data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    g0covs <- model.matrix(recharge$g0,m$data[aInd,])
    nbG0covs <- ncol(g0covs)-1
    recovs <- model.matrix(recharge$theta,m$data)
    nbRecovs <- ncol(recovs)-1
    m$data$recharge<-rep(0,nrow(m$data))
    for(i in 1:nbAnimals){
      idInd <- which(m$data$ID==unique(m$data$ID)[i])
      if(nbRecovs){
        g0 <- m$mle$g0 %*% t(g0covs[i,,drop=FALSE])
        theta <- m$mle$theta
        m$data$recharge[idInd] <- cumsum(c(g0,theta%*%t(recovs[idInd[-length(idInd)],])))
      }
    }
    if(is.null(tempCovs$recharge)) tempCovs$recharge <- mean(m$data$recharge)
    newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
    nbCovs <- nbCovs + 1
  } else {
    nbG0covs <- 0
    nbRecovs <- 0
  }

  # inverse of Hessian
  if(!is.null(m$mod$hessian)) Sigma <- m$mod$Sigma
  else Sigma <- NULL
  
  ncmean <- get_ncmean(distnames,m$conditions$fullDM,m$conditions$circularAngleMean,nbStates)
  nc <- ncmean$nc
  meanind <- ncmean$meanind
  
  tmPar <- lapply(m$mle[distnames],function(x) c(t(x)))
  parCount<- lapply(m$conditions$fullDM,ncol)
  for(i in distnames[!unlist(lapply(m$conditions$circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
  }
  parindex <- c(0,cumsum(unlist(parCount))[-length(m$conditions$fullDM)])
  names(parindex) <- distnames
  
  for(i in distnames){
    if(!is.null(m$conditions$DM[[i]])){# & m$conditions$DMind[[i]]){
      tmPar[[i]] <- m$mod$estimate[parindex[[i]]+1:parCount[[i]]]
      if(!isFALSE(m$conditions$circularAngleMean[[i]])){
        names(tmPar[[i]]) <- unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]]))))
      } else names(tmPar[[i]])<-colnames(m$conditions$fullDM[[i]])
    } else{
      if(dist[[i]] %in% angledists)
        if(!m$conditions$estAngleMean[[i]])
          tmPar[[i]] <- tmPar[[i]][-(1:nbStates)]
    }
  }
  
  Par <- list()
  lower<-list()
  upper<-list()
  se<-list()
  
  inputs <- checkInputs(nbStates,dist,tmPar,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  p<-inputs$p
  splineInputs<-getSplineDM(distnames,inputs$DM,m,tempCovs)
  covs<-splineInputs$covs
  DMinputs<-getDM(covs,splineInputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,tmPar,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
  fullDM <- DMinputs$fullDM
  #DMind <- DMinputs$DMind
  #DMinputs<-getDM(tempCovs,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,tmPar,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
  #fullDM<-DMinputs$fullDM
  
  for(i in distnames){
    tmpParNames <- p$parNames[[i]]
    tmpParNames[which(p$parNames[[i]]=="kappa")] <- "concentration"
    if(!m$conditions$DMind[[i]]){
      par <- c(w2n(m$mod$estimate,p$bounds,p$parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean[i],inputs$consensus[i],m$conditions$stationary,m$conditions$cons,fullDM,m$conditions$DMind,m$conditions$workcons,1,inputs$dist[i],m$conditions$Bndind,nc,meanind,m$covsDelta,m$conditions$workBounds,m$covsPi)[[i]])
    } else {
      par <- as.vector(t(m$mle[[i]]))
    }
    if(!(inputs$dist[[i]] %in% angledists) | (inputs$dist[[i]] %in% angledists & m$conditions$estAngleMean[[i]] & !m$conditions$Bndind[[i]])) {
      Par[[i]] <- get_CI(m$mod$estimate,par,m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],Sigma,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],m$conditions$workBounds[[i]])
    } else {
      if(!m$conditions$estAngleMean[[i]])
        Par[[i]] <- get_CI(m$mod$estimate,par[-c(1:nbStates)],m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],Sigma,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],m$conditions$workBounds[[i]])
      else {
        if(m$conditions$Bndind[[i]]){
          Par[[i]] <- CI_angle(m$mod$estimate,par,m,parindex[[i]]+1:parCount[[i]],fullDM[[i]],DMind[[i]],p$bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],Sigma,m$conditions$circularAngleMean[[i]],inputs$consensus[[i]],nbStates,alpha,tmpParNames,m$stateNames,nc[[i]],meanind[[i]],m$conditions$workBounds[[i]])
        }
      }
    }
  }
  
  mixtures <- m$conditions$mixtures

  if(nbStates>1) {
    
    # identify parameters of interest
    i2 <- tail(cumsum(unlist(parCount)),1)+1
    i3 <- i2+nbStates*(nbStates-1)*(nbCovs+1)*mixtures-1
    
    quantSup <- qnorm(1-(1-alpha)/2)
    tmpSplineInputs<-getSplineFormula(newformula,m$data,tempCovs)
    tempCovMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
    
    est<-lower<-upper<-se<-matrix(NA,nbStates*mixtures,nbStates)
    for(mix in 1:mixtures){
    
      if(is.null(recharge)){
        wpar <- m$mod$estimate[i2:i3][unique(c(m$conditions$betaCons))]
        est[(mix-1)*nbStates+1:nbStates,] <- get_gamma(wpar,tempCovMat,nbStates,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=m$conditions$workBounds$beta,mixture=mix)
        tmpSig <- Sigma[(i2:i3)[unique(c(m$conditions$betaCons))],(i2:i3)[unique(c(m$conditions$betaCons))]]
      } else {
        wpar <- c(m$mod$estimate[i2:i3][unique(c(m$conditions$betaCons))],m$mod$estimate[length(m$mod$estimate)-nbRecovs:0])
        est[(mix-1)*nbStates+1:nbStates,] <- get_gamma_recharge(wpar,tmpSplineInputs$covs,tmpSplineInputs$formula,recharge,nbStates,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=rbind(m$conditions$workBounds$beta,m$conditions$workBounds$theta),mixture=mix)
        tmpSig <- Sigma[c((i2:i3)[unique(c(m$conditions$betaCons))],length(m$mod$estimate)-nbRecovs:0),c((i2:i3)[unique(c(m$conditions$betaCons))],length(m$mod$estimate)-nbRecovs:0)]
      }
  
      if(!is.null(Sigma)){
        for(i in 1:nbStates){
          for(j in 1:nbStates){
            if(is.null(recharge)){
              dN<-numDeriv::grad(get_gamma,wpar,covs=tempCovMat,nbStates=nbStates,i=i,j=j,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=m$conditions$workBounds$beta,mixture=mix)
            } else {
              dN<-numDeriv::grad(get_gamma_recharge,wpar,covs=tmpSplineInputs$covs,formula=tmpSplineInputs$formula,recharge=recharge,nbStates=nbStates,i=i,j=j,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=rbind(m$conditions$workBounds$beta,m$conditions$workBounds$theta),mixture=mix)
            }  
            se[(mix-1)*nbStates+i,j]<-suppressWarnings(sqrt(dN%*%tmpSig%*%dN))
            lower[(mix-1)*nbStates+i,j]<-1/(1+exp(-(log(est[(mix-1)*nbStates+i,j]/(1-est[(mix-1)*nbStates+i,j]))-quantSup*(1/(est[(mix-1)*nbStates+i,j]-est[(mix-1)*nbStates+i,j]^2))*se[(mix-1)*nbStates+i,j])))#est[i,j]-quantSup*se[i,j]
            upper[(mix-1)*nbStates+i,j]<-1/(1+exp(-(log(est[(mix-1)*nbStates+i,j]/(1-est[(mix-1)*nbStates+i,j]))+quantSup*(1/(est[(mix-1)*nbStates+i,j]-est[(mix-1)*nbStates+i,j]^2))*se[(mix-1)*nbStates+i,j])))#m$mle$gamma[i,j]+quantSup*se[i,j]
          }
        }
      }
    }
    Par$gamma <- list(est=est,se=se,lower=lower,upper=upper)
    dimnames(Par$gamma$est) <- dimnames(Par$gamma$se) <- dimnames(Par$gamma$lower) <- dimnames(Par$gamma$upper) <- list(rep(m$stateNames,mixtures),m$stateNames)
    if(mixtures>1) dimnames(Par$gamma$est) <- dimnames(Par$gamma$se) <- dimnames(Par$gamma$lower) <- dimnames(Par$gamma$upper) <- list(paste0(rep(m$stateNames,mixtures),"_mix",rep(1:mixtures,each=nbStates)),m$stateNames)

    # pi
    if(mixtures>1){
      wpar<-m$mod$estimate
      pie <- matrix(wpar[i3+1:(ncol(m$covsPi)*(mixtures-1))],nrow=ncol(m$covsPi),ncol=mixtures-1)
      tmpSig <- Sigma[i3+1:(ncol(m$covsPi)*(mixtures-1)),i3+1:(ncol(m$covsPi)*(mixtures-1))]
      quantSup <- qnorm(1-(1-alpha)/2)
      lower<-upper<-se<-matrix(NA,nrow=nrow(m$covsPi),ncol=mixtures)
      est<-matrix(m$mle$pi,nrow=nrow(m$covsPi),ncol=mixtures)
      if(!is.null(Sigma)){
        for(j in 1:nrow(m$covsPi)){
          for(i in 1:mixtures){
            dN<-numDeriv::grad(get_delta,pie,covsDelta=m$covsPi[j,,drop=FALSE],i=i,workBounds=m$conditions$workBounds$pi)
            se[j,i]<-suppressWarnings(sqrt(dN%*%tmpSig%*%dN))
            lower[j,i]<-1/(1+exp(-(log(est[j,i]/(1-est[j,i]))-quantSup*(1/(est[j,i]-est[j,i]^2))*se[j,i])))
            upper[j,i]<-1/(1+exp(-(log(est[j,i]/(1-est[j,i]))+quantSup*(1/(est[j,i]-est[j,i]^2))*se[j,i])))
          }
        }
      }

      Par$pi <- list(est=est,se=se,lower=lower,upper=upper)  
      colnames(Par$pi$est) <- colnames(Par$pi$se) <- colnames(Par$pi$lower) <- colnames(Par$pi$upper) <- paste0("mix",1:mixtures)
      rownames(Par$pi$est) <- rownames(Par$pi$se) <- rownames(Par$pi$lower) <- rownames(Par$pi$upper) <- paste0("ID:",unique(m$data$ID))
    }
    
    # delta
    if(!m$conditions$stationary){
      wpar<-m$mod$estimate
      nbCovsDelta <- ncol(m$covsDelta)-1
      foo <- length(wpar)-ifelse(nbRecovs,(nbRecovs+1)+(nbG0covs+1),0)-(nbCovsDelta+1)*(nbStates-1)*mixtures
      delta <- matrix(wpar[foo+1:((nbCovsDelta+1)*(nbStates-1)*mixtures)],nrow=(nbCovsDelta+1)*mixtures,ncol=nbStates-1)
      tmpSig <- Sigma[foo+1:((nbCovsDelta+1)*(nbStates-1)*mixtures),foo+1:((nbCovsDelta+1)*(nbStates-1)*mixtures)]
      quantSup <- qnorm(1-(1-alpha)/2)
      lower<-upper<-se<-matrix(NA,nrow=nrow(m$covsDelta)*mixtures,ncol=nbStates)
      est<-matrix(m$mle$delta,nrow=nrow(m$covsDelta)*mixtures,ncol=nbStates)
      for(mix in 1:mixtures){
        if(!is.null(Sigma)){
          for(j in 1:nrow(m$covsDelta)){
            for(i in 1:nbStates){
              dN<-numDeriv::grad(get_delta,delta,covsDelta=m$covsDelta[j,,drop=FALSE],i=i,workBounds=m$conditions$workBounds$delta,mixture=mix)
              se[(mix-1)*nrow(m$covsDelta)+j,i]<-suppressWarnings(sqrt(dN%*%tmpSig%*%dN))
              lower[(mix-1)*nrow(m$covsDelta)+j,i]<-1/(1+exp(-(log(est[(mix-1)*nrow(m$covsDelta)+j,i]/(1-est[(mix-1)*nrow(m$covsDelta)+j,i]))-quantSup*(1/(est[(mix-1)*nrow(m$covsDelta)+j,i]-est[(mix-1)*nrow(m$covsDelta)+j,i]^2))*se[(mix-1)*nrow(m$covsDelta)+j,i])))#est[(mix-1)*nrow(m$covsDelta)+j,i]-quantSup*se[i]
              upper[(mix-1)*nrow(m$covsDelta)+j,i]<-1/(1+exp(-(log(est[(mix-1)*nrow(m$covsDelta)+j,i]/(1-est[(mix-1)*nrow(m$covsDelta)+j,i]))+quantSup*(1/(est[(mix-1)*nrow(m$covsDelta)+j,i]-est[(mix-1)*nrow(m$covsDelta)+j,i]^2))*se[(mix-1)*nrow(m$covsDelta)+j,i])))#est[j,i]+quantSup*se[i]
            }
          }
        }
      }
    } else {
      covs<-tempCovMat
      statFun<-function(beta,nbStates,covs,i,mixture=1){
        gamma <- trMatrix_rcpp(nbStates,beta[(mixture-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],covs,m$conditions$betaRef)[,,1]
        tryCatch(solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))[i],error = function(e) {
          "A problem occurred in the calculation of the stationary distribution."})
      }
      est <- lower <- upper <- se <- matrix(NA,nbAnimals*mixtures,nbStates)
      wpar <- m$mod$estimate[i2:i3][unique(c(m$conditions$betaCons))]
      tmpSig <- Sigma[(i2:i3)[unique(c(m$conditions$betaCons))],(i2:i3)[unique(c(m$conditions$betaCons))]]
      
      for(mix in 1:mixtures){
        delta <- statFun(matrix(wpar,nrow=(nbCovs+1)*mixtures),nbStates,covs,1:nbStates,mixture=mix)
        est[nbAnimals*(mix-1)+1:nbAnimals,] <- matrix(delta,nrow=nbAnimals,ncol=nbStates,byrow=TRUE)
        for(k in 1:nbStates){
          dN<-numDeriv::grad(statFun,matrix(wpar,nrow=(nbCovs+1)*mixtures),nbStates=nbStates,covs=covs,i=k,mixture=mix)
          se[nbAnimals*(mix-1)+1:nbAnimals,k]<-suppressWarnings(sqrt(dN%*%tmpSig%*%dN))
          lower[nbAnimals*(mix-1)+1:nbAnimals,k] <- probCI(est[nbAnimals*(mix-1)+1:nbAnimals,k],se[nbAnimals*(mix-1)+1:nbAnimals,k],quantSup,bound="lower")
          upper[nbAnimals*(mix-1)+1:nbAnimals,k] <- probCI(est[nbAnimals*(mix-1)+1:nbAnimals,k],se[nbAnimals*(mix-1)+1:nbAnimals,k],quantSup,bound="upper")
        }
      }      
    }
    Par$delta <- list(est=est,se=se,lower=lower,upper=upper)  
    
  } else {
    Par$delta <- list(est=matrix(1,nrow(m$covsDelta)),se=matrix(NA,nrow(m$covsDelta)),lower=matrix(NA,nrow(m$covsDelta)),upper=matrix(NA,nrow(m$covsDelta)))
  }
  colnames(Par$delta$est) <- colnames(Par$delta$se) <- colnames(Par$delta$lower) <- colnames(Par$delta$upper) <- m$stateNames
  rownames(Par$delta$est) <- rownames(Par$delta$se) <- rownames(Par$delta$lower) <- rownames(Par$delta$upper) <- paste0("ID:",rep(unique(m$data$ID),mixtures))
  if(mixtures>1) rownames(Par$delta$est) <- rownames(Par$delta$se) <- rownames(Par$delta$lower) <- rownames(Par$delta$upper) <- paste0("ID:",rep(unique(m$data$ID),mixtures),"_mix",rep(1:mixtures,each=nbAnimals))
  return(Par)
}

get_gamma <- function(beta,covs,nbStates,i,j,betaRef,betaCons,workBounds=NULL,mixture=1){
  tmpBeta <- rep(NA,length(betaCons))
  tmpBeta[unique(c(betaCons))] <- beta
  beta <- w2wn(matrix(tmpBeta[betaCons],nrow(betaCons),ncol(betaCons)),workBounds)
  gamma <- trMatrix_rcpp(nbStates,beta[(mixture-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],covs,betaRef)[,,1]
  gamma[i,j]
}

get_recharge <- function(g0theta,g0covs,recovs,workBounds=NULL,k=0){
  
  g0 <- w2wn(g0theta[1:ncol(g0covs)],workBounds$g0)
  theta <- w2wn(g0theta[-(1:ncol(g0covs))],workBounds$theta)
  
  g <- g0 %*% t(g0covs)
  recharge <- cumsum(c(g,theta%*%t(recovs)))
  if(k) recharge <- recharge[k]
  return(recharge)
}

get_gamma_recharge <- function(beta,covs,formula,recharge,nbStates,i,j,betaRef,betaCons,workBounds=NULL,mixture=1){
  
  recovs <- model.matrix(recharge$theta,covs)
  
  beta <- w2wn(c(beta[betaCons],beta[length(beta)-(ncol(recovs)-1):0]),workBounds)

  #g0 <- beta[length(beta)-ncol(recovs)-(ncol(g0covs)-1):0] %*% t(g0covs)
  theta <- beta[length(beta)-(ncol(recovs)-1):0]
  covs[,"recharge"] <- covs[,"recharge"] + theta%*%t(recovs) # g0  + theta%*%t(recovs)
  #covs <- matrix(covs,1)

  newcovs <- model.matrix(formula,covs)
  beta <- matrix(beta[1:(length(beta)-(ncol(recovs)))],ncol=nbStates*(nbStates-1))
  gamma <- trMatrix_rcpp(nbStates,beta[(mixture-1)*ncol(newcovs)+1:ncol(newcovs),,drop=FALSE],newcovs,betaRef)[,,1]
  gamma[i,j]
}

get_delta <- function(delta,covsDelta,i,workBounds=matrix(c(-Inf,Inf),length(delta),2,byrow=TRUE),mixture=1){
  delta <- w2wn(delta,workBounds)
  nbCovsDelta <- ncol(covsDelta)-1
  delta <- delta[(mixture-1)*(nbCovsDelta+1)+1:(nbCovsDelta+1),]
  delta <- c(rep(0,nbCovsDelta+1),delta)
  deltaXB <- covsDelta%*%matrix(delta,nrow=nbCovsDelta+1)
  expdelta <- exp(deltaXB)
  if(!is.finite(sum(expdelta))){
    tmp <- exp(Brobdingnag::as.brob(as.vector(deltaXB)))
    delta <- as.numeric(tmp/Brobdingnag::sum(tmp))
  } else {
    delta <- expdelta/sum(expdelta)
  }
  delta[i]
}

parm_list<-function(est,se,lower,upper,rnames,cnames){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  rownames(Par$est) <- rnames
  rownames(Par$se) <- rnames
  rownames(Par$lower) <- rnames
  rownames(Par$upper) <- rnames
  colnames(Par$est) <- cnames
  colnames(Par$se) <- cnames
  colnames(Par$lower) <- cnames
  colnames(Par$upper) <- cnames
  Par
}

get_CI<-function(wpar,Par,m,ind,DM,DMind,Bounds,cons,workcons,Sigma,circularAngleMean,consensus,nbStates,alpha,rnames,cnames,nc,meanind,workBounds){

  w<-wpar[ind]
  lower<-upper<-se<-rep(NA,nrow(DM))
  if(!is.null(Sigma)){
    for(k in 1:nrow(DM)){
      dN<-tryCatch(numDeriv::grad(w2nDM,w,bounds=Bounds,DM=DM,DMind=DMind,cons=cons,workcons=workcons,nbObs=1,circularAngleMean=circularAngleMean,consensus=consensus,nbStates=nbStates,k=k,nc=nc,meanind=meanind,workBounds=workBounds),error=function(e) e)
      if(inherits(dN,"error")){
        warning(dN)
        dN <- rep(NaN,length(w))
      }
      se[k]<-suppressWarnings(sqrt(dN%*%Sigma[ind,ind]%*%dN))
      lower[k] <- Par[k] - qnorm(1-(1-alpha)/2) * se[k]
      upper[k] <- Par[k] + qnorm(1-(1-alpha)/2) * se[k]
      #cn<-exp(qnorm(1-(1-alpha)/2)*sqrt(log(1+(se[k]/Par[k])^2)))
      #lower[k]<-Par[k]/cn
      #upper[k]<-Par[k]*cn
    }
  }
  est<-matrix(Par,ncol=nbStates,byrow=TRUE)
  l<-matrix(lower,ncol=nbStates,byrow=TRUE)
  u<-matrix(upper,ncol=nbStates,byrow=TRUE)  
  s<-matrix(se,ncol=nbStates,byrow=TRUE)
  out <- parm_list(est,s,l,u,rnames,cnames)
  out
}

CI_angle<-function(wpar,Par,m,ind,DM,DMind,Bounds,cons,workcons,Sigma,circularAngleMean,consensus,nbStates,alpha,rnames,cnames,nc,meanind,workBounds){
  
  w<-wpar[ind]
  lower<-upper<-se<-rep(NA,nrow(DM))
  if(!is.null(Sigma)){
    for(k in 1:nrow(DM)){
      dN<-numDeriv::grad(w2nDMangle,w,bounds=Bounds,DM=DM,DMind=DMind,cons=cons,workcons=workcons,nbObs=1,circularAngleMean=circularAngleMean,consensus=consensus,nbStates=nbStates,k=k,nc=nc,meanind=meanind,workBounds=workBounds)
      se[k]<-suppressWarnings(sqrt(dN%*%Sigma[ind,ind]%*%dN))
      lower[k] <- Par[k] - qnorm(1-(1-alpha)/2) * se[k]
      upper[k] <- Par[k] + qnorm(1-(1-alpha)/2) * se[k]
      #cn<-exp(qnorm(1-(1-alpha)/2)*sqrt(log(1+(se[k]/Par[k])^2)))
      #lower[k]<-Par[k]/cn
      #upper[k]<-Par[k]*cn
    }
  }
  est<-matrix(Par,ncol=nbStates,byrow=TRUE)
  l<-matrix(lower,ncol=nbStates,byrow=TRUE)
  u<-matrix(upper,ncol=nbStates,byrow=TRUE)  
  s<-matrix(se,ncol=nbStates,byrow=TRUE)
  out <- parm_list(est,s,l,u,rnames,cnames)
  out
}
  
w2nDMangle<-function(w,bounds,DM,DMind,cons,workcons,nbObs,circularAngleMean,consensus,nbStates,k,nc,meanind,workBounds){
  
  bounds[,1] <- -Inf
  bounds[which(bounds[,2]!=1),2] <- Inf
    
  foo <- length(w) - nbStates + 1
  x <- w[(foo - nbStates):(foo - 1)]
  y <- w[foo:length(w)]
  angleMean <- Arg(x + (0+1i) * y)
  kappa <- sqrt(x^2 + y^2)
  w[(foo - nbStates):(foo - 1)] <- angleMean
  w[foo:length(w)] <- kappa
  
  w2nDM(w,bounds,DM,DMind,cons,workcons,nbObs,circularAngleMean,consensus,nbStates,k,nc,meanind,workBounds)
}
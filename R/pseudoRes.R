
#' Pseudo-residuals
#'
#' The pseudo-residuals of momentuHMM models, as described in Zucchini and McDonad (2009).
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object.
#'
#' @return A list of psuedo-residuals for each data stream (e.g., 'stepRes', 'angleRes')
#'
#' @details If some turning angles in the data are equal to pi, the corresponding pseudo-residuals
#' will not be included. Indeed, given that the turning angles are defined on (-pi,pi], an angle of pi
#' results in a pseudo-residual of +Inf (check Section 6.2 of reference for more information on the
#' computation of pseudo-residuals).
#' 
#' Note that pseudo-residuals for multiple imputation analyses are based on pooled parameter 
#' estimates and the means of the data values across all imputations.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' res <- pseudoRes(m)
#' qqnorm(res$stepRes)
#' qqnorm(res$angleRes)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export
#' @importFrom stats integrate qnorm

pseudoRes <- function(m)
{
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(m)) m <- m$miSum

  data <- m$data
  nbObs <- nrow(data)
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  
  if(is.miSum(m)){
    warning('pseudo-residuals are based on pooled parameter estimates and mean data values across multiple imputations...')
    Par <- lapply(m$Par$real,function(x) x$est)
    for(i in distnames){
      if(!is.null(m$conditions$DM[[i]]))
        Par[[i]] <- m$Par$beta[[i]]$est
      else if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]])
        Par[[i]] <- Par[[i]][-1,]
    }
    Par<-lapply(Par,function(x) c(t(x)))
    Par<-Par[distnames]
    beta <- m$Par$beta$beta$est
    delta <- m$Par$real$delta$est
    inputs <- checkInputs(nbStates,m$conditions$dist,Par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
    p <- inputs$p
    DMinputs<-getDM(data,inputs$DM,m$conditions$dist,nbStates,p$parNames,p$bounds,Par,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$circularAngleMean)
    m$conditions$fullDM <- DMinputs$fullDM
    m$mod$estimate <- n2w(Par,p$bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
  }
  
  Fun <- lapply(dist,function(x) paste("p",x,sep=""))
  for(j in which(dist %in% angledists)){
    Fun[[j]] <- paste0("d",dist[[j]])
    if(length(which(data[[distnames[j]]]==pi))>0)
      message("Note: Some ",distnames[j],"s are equal to pi, and the corresponding pseudo-residuals are not included")
  }

  # forward log-probabilities
  la <- logAlpha(m)
  
  # identify covariates
  formula<-m$conditions$formula
  stateForms<- terms(formula, specials = paste0("state",1:nbStates))
  newformula<-formula
  if(nbStates>1){
    if(length(unlist(attr(stateForms,"specials")))){
      newForm<-attr(stateForms,"term.labels")[-unlist(attr(stateForms,"specials"))]
      for(i in 1:nbStates){
        if(!is.null(attr(stateForms,"specials")[[paste0("state",i)]])){
          for(j in 1:(nbStates-1)){
            newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("state",i)]]]))
          }
        }
      }
      newformula<-as.formula(paste("~",paste(newForm,collapse="+")))
    }
    formulaStates<-stateFormulas(newformula,nbStates*(nbStates-1),spec="betaCol")
    if(length(unlist(attr(terms(newformula, specials = c(paste0("betaCol",1:(nbStates*(nbStates-1))),"cosinor")),"specials")))){
      allTerms<-unlist(lapply(formulaStates,function(x) attr(terms(x),"term.labels")))
      newformula<-as.formula(paste("~",paste(allTerms,collapse="+")))
      formterms<-attr(terms.formula(newformula),"term.labels")
    } else {
      formterms<-attr(terms.formula(newformula),"term.labels")
      newformula<-formula
    }
  }
  covs <- model.matrix(newformula,data)
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
  }
  
  par <- w2n(m$mod$estimate,m$conditions$bounds,lapply(m$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,nbObs,dist,m$conditions$Bndind,nc,meanind)
  
  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  genRes <- list()
  for(j in distnames){
    genRes[[paste0(j,"Res")]] <- rep(NA,nbObs)
    pgenMat <- matrix(NA,nbObs,nbStates)
    sp <- par[[j]]
    genInd <- which(!is.na(data[[j]]))
    zeroInflation <- m$conditions$zeroInflation[[j]]
    oneInflation <- m$conditions$oneInflation[[j]]
  
    for(state in 1:nbStates) {
      
      genPar <- sp
      
      if(!(dist[[j]] %in% angledists)){
        
        genArgs <- list(data[[j]][genInd])
        
        zeromass <- 0
        onemass <- 0
        if(zeroInflation | oneInflation) {
          if(zeroInflation) zeromass <- genPar[nrow(genPar)-nbStates*oneInflation-nbStates+state,genInd]
          if(oneInflation) onemass <- genPar[nrow(genPar)-nbStates+state,genInd]
          genPar <- genPar[-(nrow(genPar)-(nbStates*(zeroInflation+oneInflation)-1):0),]
        }
        for(k in 1:(nrow(genPar)/nbStates))
          genArgs[[k+1]] <- genPar[(k-1)*nbStates+state,genInd]
        
        if(dist[[j]]=="gamma") {
          shape <- genArgs[[2]]^2/genArgs[[3]]^2
          scale <- genArgs[[3]]^2/genArgs[[2]]
          genArgs[[2]] <- shape
          genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        
        if(zeroInflation | oneInflation) {
          if(zeroInflation & !oneInflation){
            pgenMat[genInd,state] <- ifelse(data[[j]][genInd]==0,
                                      zeromass, # if gen==0
                                      (1-zeromass)*do.call(Fun[[j]],genArgs)) # if gen != 0
          } else if(oneInflation & !zeroInflation){
            pgenMat[genInd,state] <- ifelse(data[[j]][genInd]==1,
                                      onemass, # if gen==0
                                      (1-onemass)*do.call(Fun[[j]],genArgs)) # if gen != 1          
          } else {
            pgenMat[genInd,state][data[[j]][genInd]==0] <- zeromass[data[[j]][genInd]==0]
            pgenMat[genInd,state][data[[j]][genInd]==1] <- (1.-zeromass[data[[j]][genInd]==1]) * onemass[data[[j]][genInd]==1]
            pgenMat[genInd,state][data[[j]][genInd]>0 & data[[j]][genInd]<1] <- (1.-zeromass[data[[j]][genInd]>0 & data[[j]][genInd]<1]) * (1.-onemass[data[[j]][genInd]>0 & data[[j]][genInd]<1]) * do.call(Fun[[j]],genArgs)[data[[j]][genInd]>0 & data[[j]][genInd]<1] # if gen !=0 and gen!=1
          }
        }
        else pgenMat[genInd,state] <- do.call(Fun[[j]],genArgs)
        
        
        #pgenMat[genInd,state] <- zeromass+(1-zeromass)*do.call(Fun[[j]],genArgs)
        #for(i in 1:nbObs) {
        #  if(!is.na(data[[j]][i])) {
        #    genArgs[[1]] <- data[[j]][i]
        #    pgenMat[i,state] <- zeromass+(1-zeromass)*do.call(Fun[[j]],genArgs)
        #  }
        #}
      } else {
        
        genpiInd <- which(data[[j]]!=pi & !is.na(data[[j]]))
        
        genArgs <- list(Fun[[j]],-pi,data[[j]][1]) # to pass to function "integrate" below
  
        for(i in genpiInd){
          genArgs[[3]]<-data[[j]][i]
          for(k in 1:(nrow(genPar)/nbStates))
            genArgs[[k+3]] <- genPar[(k-1)*nbStates+state,i]
          
          pgenMat[i,state] <- do.call(integrate,genArgs)$value
        }
        #for(i in 1:nbObs) {
        #  if(!is.na(data[[j]][i])) {
        #    # angle==pi => residual=Inf
        #    if(data[[j]][i]!=pi) {
        #      genArgs[[3]] <- data[[j]][i]
        #      pgenMat[i,state] <- do.call(integrate,genArgs)$value
        #    }
        #  }
        #}
      }
    }
  
    if(!is.na(data[[j]][1]))
      genRes[[paste0(j,"Res")]][1] <- qnorm((delta%*%trMat[,,1])%*%pgenMat[1,])

    for(i in 2:nbObs) {
      gamma <- trMat[,,i]
      c <- max(la[i-1,]) # cancels below ; prevents numerical errors
      a <- exp(la[i-1,]-c)
  
      if(!is.na(data[[j]][i]))
        genRes[[paste0(j,"Res")]][i] <-qnorm(t(a)%*%(gamma/sum(a))%*%pgenMat[i,])
    }
  }

  return(genRes)
}

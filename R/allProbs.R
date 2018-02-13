
#' Matrix of all probabilities
#'
#' Used in functions \code{\link{viterbi}}, \code{\link{logAlpha}}, \code{\link{logBeta}}.
#'
#' @param m Object \code{\link{momentuHMM}} or \code{\link{miSum}}.
#' @param nbStates Number of states of the HMM.
#'
#' @return Matrix of all probabilities.
#'
#' @examples
#' \dontrun{
#' P <- momentuHMM:::allProbs(m=example$m,nbStates=2)
#' }
#' @importFrom LaplacesDemon dbern

allProbs <- function(m,nbStates)
{
  
  if(!is.momentuHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM) or miSum object (as returned by MIpool)")
  
  m <- delta_bc(m)
  
  data <- m$data
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  zeroInflation <- m$conditions$zeroInflation
  oneInflation <- m$conditions$oneInflation
  nbObs <- nrow(data)
  
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
  nbCovs <- ncol(model.matrix(newformula,data))-1 # substract intercept column
  
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) {
      meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  
  consensus <- vector('list',length(distnames))
  names(consensus) <- distnames
  for(i in distnames){
    consensus[[i]] <- (dist[[i]]=="vmConsensus")
  }
  dist <- lapply(dist,function(x) gsub("Consensus","",x))

  par <- w2n(m$mod$estimate,m$conditions$bounds,lapply(m$conditions$fullDM,function(x) nrow(x)/nbStates),nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,consensus,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,nbObs,dist,m$conditions$Bndind,nc,meanind,m$covsDelta,m$conditions$workBounds)
  
  Fun <- lapply(dist,function(x) paste("d",x,sep=""))
  

  probs <- matrix(1,nrow=nbObs,ncol=nbStates)
  
  for(i in distnames){
  
    genInd <- which(!is.na(data[[i]]))
    sp <- par[[i]]
  
    for(state in 1:nbStates) {
      genPar <- sp
      genProb <- rep(1,nbObs)
      genFun <- Fun[[i]]
      
      # Constitute the lists of state-dependent parameters for the step and angle
      genArgs <- list(data[[i]][genInd])
      
      zeromass <- 0
      onemass <- 0
      if(zeroInflation[[i]] | oneInflation[[i]]) {
        if(zeroInflation[[i]]) zeromass <- genPar[nrow(genPar)-nbStates*oneInflation[[i]]-nbStates+state,genInd]
        if(oneInflation[[i]]) onemass <- genPar[nrow(genPar)-nbStates+state,genInd]
        genPar <- genPar[-(nrow(genPar)-(nbStates*(zeroInflation[[i]]+oneInflation[[i]])-1):0),]
      }
  
      for(j in 1:(nrow(genPar)/nbStates))
        genArgs[[j+1]] <- genPar[(j-1)*nbStates+state,genInd]
  
      # conversion between mean/sd and shape/scale if necessary
      if(dist[[i]]=="gamma") {
        shape <- genArgs[[2]]^2/genArgs[[3]]^2
        scale <- genArgs[[3]]^2/genArgs[[2]]
        genArgs[[2]] <- shape
        genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
      }
      if(zeroInflation[[i]] | oneInflation[[i]]) {
        if(zeroInflation[[i]] & !oneInflation[[i]]){
          genProb[genInd] <- ifelse(data[[i]][genInd]==0,
                                      zeromass, # if gen==0
                                      (1-zeromass)*do.call(genFun,genArgs)) # if gen != 0
        } else if(oneInflation[[i]] & !zeroInflation[[i]]){
          genProb[genInd] <- ifelse(data[[i]][genInd]==1,
                                    onemass, # if gen==0
                                    (1-onemass)*do.call(genFun,genArgs)) # if gen != 1          
        } else {
          genProb[genInd][data[[i]][genInd]==0] <- zeromass[data[[i]][genInd]==0]
          genProb[genInd][data[[i]][genInd]==1] <- onemass[data[[i]][genInd]==1]
          genProb[genInd][data[[i]][genInd]>0 & data[[i]][genInd]<1] <- (1.-zeromass[data[[i]][genInd]>0 & data[[i]][genInd]<1]-onemass[data[[i]][genInd]>0 & data[[i]][genInd]<1]) * do.call(genFun,genArgs)[data[[i]][genInd]>0 & data[[i]][genInd]<1] # if gen !=0 and gen!=1
        }
      }
      else genProb[genInd] <- do.call(genFun,genArgs)
  
      probs[,state] <- probs[,state]*genProb;
    }
  }
  if(!is.null(m$knownStates)) {
    for (i in which(!is.na(m$knownStates))) {
      prob <- probs[i, m$knownStates[i]]
      probs[i, ] <- 0
      probs[i, m$knownStates[i]] <- prob
    }
  }
  return(probs)
}

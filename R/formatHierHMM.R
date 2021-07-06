
#' Convert hierarchical HMM structure to a conventional HMM
#' 
#' @param data \code{\link{momentuHierHMMData}} object or a data frame containing the data streams and covariates. 
#' @param hierStates A hierarchical data structure \code{\link[data.tree]{Node}} for the states ('state').  See \code{\link{fitHMM}}. 
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). See \code{\link{fitHMM}}. 
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the transition probabilities at each level of the hierarchy ('beta'). See \code{\link{fitHMM}}. 
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the initial distribution at each level of the hierarchy ('delta'). See \code{\link{fitHMM}}. 
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). See \code{\link{fitHMM}}. Default: \code{NULL} (only hierarchical-level effects, with no covariate effects). 
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). See \code{\link{fitHMM}}. Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). See \code{\link{fitHMM}}. Default: \code{mixtures=1}.  
#' @param workBounds A list with elements named \code{'beta'} and/or \code{'delta'}, where each element is a hierarchical data structure \code{\link[data.tree]{Node}} indicating t.p.m. and initial distribution working parameter bounds ('workBounds') for parameters in \code{hierBeta} and \code{hierDelta}, respectively.
#' @param betaCons A hierarchical data structure \code{\link[data.tree]{Node}} indicating t.p.m. constraints ('betaCons') among parameters in \code{hierBeta} at each level of the hierarchy.
#' @param deltaCons A hierarchical data structure \code{\link[data.tree]{Node}} indicating initial distribution constraints ('deltaCons') among parameters in \code{hierDelta} at each level of the hierarchy.
#' @param fixPar A list with elements named \code{'beta'} and/or \code{'delta'}, where each element is a hierarchical data structure \code{\link[data.tree]{Node}} indicating t.p.m. and initial distribution parameters in \code{hierBeta} and \code{hierDelta}, respectively, which are assumed known.
#' @param checkData logical indicating whether or not to check the suitability of \code{data} for the specified hierarchy. Ignored unless \code{data} is provided. Default: TRUE.
#' 
#' @return A list of arguments needed for specifying a hierarchical HMM as a conventional HMM in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}, including:
#' \item{nbStates}{See \code{\link{fitHMM}}.}
#' \item{dist}{See \code{\link{fitHMM}}.}
#' \item{formula}{See \code{\link{fitHMM}}.}
#' \item{formulaDelta}{See \code{\link{fitHMM}}.}
#' \item{beta0}{See \code{\link{fitHMM}}.}
#' \item{delta0}{See \code{\link{fitHMM}}.}
#' \item{betaRef}{See \code{\link{fitHMM}}.}
#' \item{betaCons}{See \code{\link{fitHMM}}.}
#' \item{deltaCons}{See \code{\link{fitHMM}}.}
#' \item{fixPar}{See \code{\link{fitHMM}}.}
#' \item{workBounds}{See \code{\link{fitHMM}}.}
#' \item{stateNames}{See \code{\link{fitHMM}}.}
#' 
#' @export
# #' @importFrom data.tree Node Get Do ToDataFrameTypeCol Traverse Aggregate AreNamesUnique isRoot isLeaf Clone
#' @importFrom stats terms
formatHierHMM <- function(data,hierStates,hierDist,
                          hierBeta=NULL,hierDelta=NULL,
                          hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,
                          workBounds=NULL,betaCons=NULL,deltaCons=NULL,fixPar=NULL,
                          checkData=TRUE){
  
  installDataTree()
  
  if(is.null(data)) checkData <- FALSE
  
  hStates <- nbHierStates(hierStates)
  nbStates <- hStates$nbStates
  stateNames <- hStates$stateNames

  dist <- getHierDist(hierDist,data,checkData)
  nbLevels <- 2*hierDist$count - 1 #hierStates$height #ncol(data.tree::ToDataFrameTypeCol(hierStates))

  whierFormulaDelta <- checkHierFormula(data,hierFormulaDelta,hierStates,hierDist,checkData,what="formulaDelta")
  formulaDelta <- whierFormulaDelta$level1$formulaDelta
  whierFormulaDelta$RemoveChild("level1")
  
  whierFormula <- checkHierFormula(data,hierFormula,hierStates,hierDist,checkData)
  for(j in names(whierFormulaDelta$children)){
    whierFormula$AddChild(paste0(j,"i"),formula=whierFormulaDelta[[j]]$formulaDelta)
  }
  
  newForm <- formatHierFormula(data,whierFormula,hierStates)
  formula <- newForm$formula
  #data <- newForm$data
  recharge <- newForm$recharge
  
  #if(length(recharge)>1) stop("sorry, only 1 recharge model is currently supported")
  
  # set t.p.m. reference states based on top level
  betaRef <- rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))

  beta0 <- delta0 <- NULL
  
  if(!is.null(data)){
    
    if(!is.null(recharge)){
      reParms <- getRechargeParms(recharge,data,hierBeta)
      g0 <- reParms$g0
      theta <- reParms$theta
      rePar <- list(g0=g0,theta=theta)
    } else rePar <- NULL
    
    reForm <- formatRecharge(nbStates,formula,betaRef,data=data,par=rePar)
    newformula <- reForm$newformula
    coords <- attr(data,"coords")
    coordLevel <- attr(data,"coordLevel")
    data <- cbind(data,reForm$newdata)
    attr(data,"coords") <- coords
    attr(data,"coordLevel") <- coordLevel
    aInd <- reForm$aInd
    covs <- stats::model.matrix(newformula,data)
    nbCovs <- ncol(covs)
    
    hFixPar <- list()
    hFixPar$beta <- matrix(NA,nbCovs*mixtures,nbStates*(nbStates-1))
    rownames(hFixPar$beta) <- paste0(rep(colnames(covs),mixtures),"_mix",rep(1:mixtures,each=length(colnames(covs))))
    colnames(hFixPar$beta) <- c(sapply(1:nbStates,function(x) paste(rep(x,each=nbStates-1),"->",1:nbStates)[-betaRef[x]]))
    
    hBetaCons <- matrix(1:(nbCovs*mixtures*nbStates*(nbStates-1)),nbCovs*mixtures,nbStates*(nbStates-1))
    dimnames(hBetaCons)<-dimnames(hFixPar$beta)
    
    covsDelta <- stats::model.matrix(formulaDelta,data[aInd,,drop=FALSE])
    nbCovsDelta <- ncol(covsDelta)
    hFixPar$delta <- matrix(NA,nbCovsDelta*mixtures,nbStates-1,byrow=TRUE)
    colnames(hFixPar$delta) <- paste("state",2:nbStates)
    rownames(hFixPar$delta) <- paste0(rep(colnames(covsDelta),mixtures),"_mix",rep(1:mixtures,each=length(colnames(covsDelta))))
    
    hDeltaCons <- matrix(1:(nbCovsDelta*mixtures*(nbStates-1)),nbCovsDelta*mixtures,(nbStates-1))
    dimnames(hDeltaCons)<-dimnames(hFixPar$delta)
    #if(any(betaRef==1)){
    for(mix in 1:mixtures){
      hFixPar$delta[paste0("(Intercept)","_mix",mix),(2:nbStates)[-(betaRef-1)]-1] <- -1.e+10
      hDeltaCons[paste0("(Intercept)","_mix",mix),(2:nbStates)[-(betaRef-1)]-1] <- min(hDeltaCons[paste0("(Intercept)","_mix",mix),(2:nbStates)[-(betaRef-1)]-1])
      if(nbCovsDelta>1) {
        hFixPar$delta[paste0(colnames(covsDelta)[which(colnames(covsDelta)!="(Intercept)")],"_mix",mix),(2:nbStates)[-(betaRef-1)]-1] <- 0
        hDeltaCons[paste0(colnames(covsDelta)[which(colnames(covsDelta)!="(Intercept)")],"_mix",mix),(2:nbStates)[-(betaRef-1)]-1] <- min(hDeltaCons[paste0(colnames(covsDelta)[which(colnames(covsDelta)=="(Intercept)")],"_mix",mix),(2:nbStates)[-(betaRef-1)]-1])
      }
    }
    #} else {
    #  if(is.null(workBounds$delta)){
    #    deltaLower <- matrix(-Inf,nbCovsDelta*mixtures,nbStates-1,dimnames=list(rownames(hFixPar$delta)))
    #    deltaUpper <- matrix( Inf,nbCovsDelta*mixtures,nbStates-1,dimnames=list(rownames(hFixPar$delta)))
    #    for(mix in 1:mixtures){
    #     deltaLower[paste0("level",levels(data$level)[1],"_mix",mix),betaRef-1] <- 100
    #     hFixPar$delta[paste0("level",levels(data$level)[1],"_mix",mix),-(betaRef-1)] <- 0
    #    }
    #    workBoundsDelta <- cbind(c(deltaLower),c(deltaUpper))
    #    if(is.null(workBounds)) workBounds <- list()
    #    workBounds$delta <- workBoundsDelta
    #  }
    #}
    
    hierStates$Do(function(node) node$levelStates <- node$Get("state",traversal="level",filterFun=data.tree::isLeaf))
    
    #top level t.p.m. (must map all states at bottom level back to top level)
    j <- 1
    t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
    names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
    intCov <- c(paste0("level",levels(data$level)[2*j-1]),paste0("I((level == \"",levels(data$level)[2*j-1],"\") * 1)"))
    intCov <- intCov[which(intCov %in% colnames(covs))]
    ointCov <- otherCovs <- list()
    for(ol in levels(data$level)[(2*j):nbLevels]){
      ointCov[[ol]] <- c(paste0("level",ol),paste0("I((level == \"",ol,"\") * 1)"))
      ointCov[[ol]] <- ointCov[[ol]][which(ointCov[[ol]] %in% colnames(covs))]
    }
    levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-1],"\")"),colnames(covs),fixed=TRUE)]
    levelCovs <- levelCovs[which(!(levelCovs %in% intCov))]
    alloCovs <- colnames(covs)[-which(colnames(covs) %in% c(intCov,levelCovs,unlist(ointCov)))]
    fInd <- which(unlist(lapply(data,function(x) inherits(x,"factor"))))
    factorCovs <- names(data)[fInd]
    factorCovs <- factorCovs[-which(factorCovs=="level")]
    nfInd <- names(data)[-fInd]
    for(ol in levels(data$level)[(2*j):nbLevels]){
      otherCovs[[ol]] <- alloCovs[grepl(paste0("I((level == \"",ol,"\")"),alloCovs,fixed=TRUE)]
      oInt <- find_intercept(ointCov[[ol]],otherCovs[[ol]],factorCovs,nfInd,data)
      ointCov[[ol]] <- oInt$intCov
      otherCovs[[ol]] <- oInt$otherCovs
    }
    ointCov <- unlist(ointCov)
    otherCovs <- unlist(otherCovs)
    Int <- find_intercept(intCov,levelCovs,factorCovs,nfInd,data)
    intCov <- Int$intCov
    levelCovs <- Int$otherCovs

    for(k in names(t)){
      tt <- data.tree::Traverse(t[[k]],filterFun=function(x) x$level==j+1)
      withinConstr <- acrossConstr <- list()
      if(length(tt)){
        for(h in 1:length(tt)){
          levelStates <- tt[[h]]$Get("state",filterFun = data.tree::isLeaf)
          stateInd <- tt[[h]]$Get("levelStates",filterFun=data.tree::isLeaf)
          withinConstr[[h]] <- match(paste0(rep(levelStates,each=length(stateInd))," -> ",stateInd),colnames(hFixPar$beta),nomatch=0)
          if(any(withinConstr[[h]])){
            for(mix in 1:mixtures){
              hBetaCons[paste0(intCov,"_mix",mix),withinConstr[[h]]] <- min(hBetaCons[paste0(intCov,"_mix",mix),withinConstr[[h]]])
              for(lcovs in levelCovs){
                hBetaCons[paste0(lcovs,"_mix",mix),withinConstr[[h]]] <- min(hBetaCons[paste0(lcovs,"_mix",mix),withinConstr[[h]]])
              }
              hFixPar$beta[paste0(intCov,"_mix",mix),withinConstr[[h]]] <- -1.e+10
              for(lcovs in levelCovs){
                hFixPar$beta[paste0(lcovs,"_mix",mix),withinConstr[[h]]] <- ifelse(length(intCov),0,-1.e+10) # if there's not intercept then all should be -1.e+10
              }
            }
          }
          # constrain across transitions to reference states
          levelStates <- lapply(tt[(1:length(tt))[-h]],function(x) x$Get("state",filterFun = data.tree::isLeaf))
          for(ll in 1:length(levelStates)){
            acrossConstr[[h]] <- match(paste0(rep(levelStates[[ll]],each=length(stateInd[-1]))," -> ",stateInd[-1]),colnames(hFixPar$beta),nomatch=0)
            acrossRef <- match(paste0(levelStates[[ll]]," -> ",stateInd[1]),colnames(hFixPar$beta),nomatch=0)
            for(mix in 1:mixtures){
              if(any(acrossConstr[[h]])){
                hBetaCons[paste0(intCov,"_mix",mix),acrossConstr[[h]]] <- min(hBetaCons[paste0(intCov,"_mix",mix),acrossConstr[[h]]])
                for(lcovs in levelCovs){
                  hBetaCons[paste0(lcovs,"_mix",mix),acrossConstr[[h]]] <- min(hBetaCons[paste0(lcovs,"_mix",mix),acrossConstr[[h]]])
                }
                hBetaCons[paste0(ointCov,"_mix",mix),acrossConstr[[h]]] <- min(hBetaCons[paste0(ointCov,"_mix",mix),acrossConstr[[h]]])
                for(ocovs in otherCovs){
                  hBetaCons[paste0(ocovs,"_mix",mix),acrossConstr[[h]]] <- min(hBetaCons[paste0(ocovs,"_mix",mix),acrossConstr[[h]]])
                }
              }
              if(any(acrossRef)){
                hBetaCons[paste0(intCov,"_mix",mix),acrossRef] <- min(hBetaCons[paste0(intCov,"_mix",mix),acrossRef])
                for(lcovs in levelCovs){
                  hBetaCons[paste0(lcovs,"_mix",mix),acrossRef] <- min(hBetaCons[paste0(lcovs,"_mix",mix),acrossRef])
                }
                hBetaCons[paste0(ointCov,"_mix",mix),acrossRef] <- min(hBetaCons[paste0(ointCov,"_mix",mix),acrossRef])
                for(ocovs in otherCovs){
                  hBetaCons[paste0(ocovs,"_mix",mix),acrossRef] <- min(hBetaCons[paste0(ocovs,"_mix",mix),acrossRef])
                }
              }
              if(any(acrossConstr[[h]])){
                hFixPar$beta[paste0(intCov,"_mix",mix),acrossConstr[[h]]] <- -1.e+10
                for(lcovs in levelCovs){
                  hFixPar$beta[paste0(lcovs,"_mix",mix),acrossConstr[[h]]] <- 0
                }
                hFixPar$beta[paste0(ointCov,"_mix",mix),acrossConstr[[h]]] <- -1.e+10
                for(ocovs in otherCovs){
                  hFixPar$beta[paste0(ocovs,"_mix",mix),acrossConstr[[h]]] <- 0
                }
              }
              if(any(acrossRef)){
                hFixPar$beta[paste0(ointCov,"_mix",mix),acrossRef] <- -1.e+10
                for(ocovs in otherCovs){
                  hFixPar$beta[paste0(ocovs,"_mix",mix),acrossRef] <- 0
                }
              }          
            }
          }
        }
      }
    }
    
    levelStateFun <- function(node) {
      if (node$parent$isRoot) return (node$levelStates)
      parent <- node$parent
      return ( levelStateFun(parent))
    }
    
    betaLower <- matrix(-Inf,nbCovs*mixtures,nbStates*(nbStates-1),dimnames = dimnames(hBetaCons))
    betaUpper <- matrix( Inf,nbCovs*mixtures,nbStates*(nbStates-1),dimnames = dimnames(hBetaCons))
    
    for(j in 2:(hierStates$height-1)){
      
      t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
      names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
      
      intCov <- c(paste0("level",levels(data$level)[2*j-2]),paste0("I((level == \"",levels(data$level)[2*j-2],"\") * 1)"))
      intCov <- intCov[which(intCov %in% colnames(covs))]
      levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-2],"\")"),colnames(covs),fixed=TRUE)]
      levelCovs <- levelCovs[which(!(levelCovs %in% intCov))]
      Int <- find_intercept(intCov,levelCovs,factorCovs,nfInd,data)
      intCov <- Int$intCov
      levelCovs <- Int$otherCovs
      
      #initial distribution       
      for(k in names(t)){
        levelStates <- t[[k]]$Get("state",filterFun = data.tree::isLeaf)
        allStates <- levelStateFun(t[[k]])
        fromState <- data.tree::Aggregate(t[[k]],"state",min)
        tt <- data.tree::Traverse(t[[k]],filterFun=function(x) x$level==j+1)
        toStates <- unlist(lapply(tt,function(x) data.tree::Aggregate(x,"state",min)))#data.tree::Aggregate(tt[[h]],"state",min)
        if(is.null(toStates)) toStates <- t[[k]]$state
        allTrans <- paste0(rep(levelStates,each=length(allStates[which(!(allStates %in% betaRef[toStates]))]))," -> ",allStates[which(!(allStates %in% betaRef[toStates]))])
        initConstr <- match(allTrans[which(!(allTrans %in% paste0(rep(fromState,each=length(toStates))," -> ",toStates)))],colnames(hFixPar$beta),nomatch=0)
        refConstr <- match(allTrans[which(allTrans %in% paste0(rep(fromState,each=length(toStates))," -> ",toStates))],colnames(hFixPar$beta),nomatch=0)
        for(mix in 1:mixtures){
          if(any(initConstr)){
            hBetaCons[paste0(intCov,"_mix",mix),initConstr] <- min(hBetaCons[paste0(intCov,"_mix",mix),initConstr])
            for(lcovs in levelCovs){
              hBetaCons[paste0(lcovs,"_mix",mix),initConstr] <- min(hBetaCons[paste0(lcovs,"_mix",mix),initConstr])
            }
            hFixPar$beta[paste0(intCov,"_mix",mix),initConstr] <- -1.e+10
            for(lcovs in levelCovs){
              hFixPar$beta[paste0(lcovs,"_mix",mix),initConstr] <- 0
            }
          }
          if(!any(betaRef[fromState] %in% toStates)){
            #warning("to state = ",toStates,"; betaRef[fromState] = ",betaRef[fromState])
            if(any(refConstr)){
              if(sum(refConstr>0)>1){
                betaLower[paste0(intCov,"_mix",mix),refConstr] <- 500
                for(lcovs in levelCovs){
                  betaLower[paste0(lcovs,"_mix",mix),refConstr] <- -(500-100)/length(levelCovs) # don't let XB be less than 100
                }
              } else if(!length(levelCovs)) hFixPar$beta[paste0(intCov,"_mix",mix),refConstr] <- 500
            }
          }
        }
      }
      
      intCov <- c(paste0("level",levels(data$level)[2*j-1]),paste0("I((level == \"",levels(data$level)[2*j-1],"\") * 1)"))
      intCov <- intCov[which(intCov %in% colnames(covs))]
      levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-1],"\")"),colnames(covs),fixed=TRUE)]
      levelCovs <- levelCovs[which(!(levelCovs %in% intCov))]
      Int <- find_intercept(intCov,levelCovs,factorCovs,nfInd,data)
      intCov <- Int$intCov
      levelCovs <- Int$otherCovs

      # t.p.m.
      for(k in names(t)){  
        levelStates <- t[[k]]$Get("state",filterFun = data.tree::isLeaf)
        allStates <- levelStateFun(t[[k]])
        tt <- data.tree::Traverse(t[[k]],filterFun=function(x) x$level==j+1)
        fromStates <- toStates <- unlist(lapply(tt,function(x) data.tree::Aggregate(x,"state",min)))
        if(is.null(toStates)) toStates <- fromStates <- t[[k]]$state
        allTrans <- paste0(rep(levelStates,each=length(allStates[which(!(allStates %in% betaRef[toStates]))]))," -> ",allStates[which(!(allStates %in% betaRef[toStates]))])
        initConstr <- match(allTrans[which(!(allTrans %in% paste0(rep(fromStates,each=length(toStates))," -> ",toStates)))],colnames(hFixPar$beta),nomatch=0)
        refConstr <- match(allTrans[which(allTrans %in% paste0(rep(fromStates,each=length(toStates))," -> ",toStates))],colnames(hFixPar$beta),nomatch=0)
        for(mix in 1:mixtures){
          if(any(initConstr)){
            hBetaCons[paste0(intCov,"_mix",mix),initConstr] <- min(hBetaCons[paste0(intCov,"_mix",mix),initConstr])
            for(lcovs in levelCovs){
              hBetaCons[paste0(lcovs,"_mix",mix),initConstr] <- min(hBetaCons[paste0(lcovs,"_mix",mix),initConstr])
            }
            hFixPar$beta[paste0(intCov,"_mix",mix),initConstr] <- -1.e+10
            for(lcovs in levelCovs){
              hFixPar$beta[paste0(lcovs,"_mix",mix),initConstr] <- 0
            }
          }
          if(!any(betaRef[fromStates] %in% toStates)){
            #warning("to state = ",toStates,"; betaRef[fromStates] = ",betaRef[fromStates])
            if(any(refConstr)) {
              if(sum(refConstr>0)>1){
                betaLower[paste0(intCov,"_mix",mix),refConstr] <- 500
                for(lcovs in levelCovs){
                  betaLower[paste0(lcovs,"_mix",mix),refConstr] <- -(500-100)/length(levelCovs) # don't let XB be less than 100
                }
              } else if(!length(levelCovs)) hFixPar$beta[paste0(intCov,"_mix",mix),refConstr] <- 500
            }
          }
        }
      }
    }
    fixInd <- which(hFixPar$beta==-1.e+10)
    if(length(fixInd)) hBetaCons[fixInd] <- fixInd[1]
    
    Pi <- NULL
    if(!is.null(hierBeta)){
      if(mixtures>1){
        #if(!is.list(hierBeta) || !all(names(hierBeta) %in% c("beta","pi","g0","theta"))) stop("hierBeta must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
        if(!is.list(hierBeta) || !all(names(hierBeta) %in% c("beta","pi"))) stop("hierBeta must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
        Pi <- hierBeta[["pi"]]
      }
      if(!is.null(recharge)){
        #  if(!is.list(hierBeta) || !all(names(hierBeta) %in% c("beta","pi","g0","theta"))) stop("hierBeta must be a list with elements named 'beta', 'g0', and/or 'theta' when including a recharge model")
        #  g0 <- hierBeta$g0
        #  theta <- hierBeta$theta
        if(is.list(hierBeta) && any(names(hierBeta) %in% c("g0","theta"))) stop("for hierarchical recharge models, initial values for g0 and theta must be specified within the corresponding levels of hierBeta")
      }
      if(is.list(hierBeta)) hierBeta <- hierBeta$beta
      whierBeta <- data.tree::Clone(hierBeta)
    } else whierBeta <- NULL
    
    if(!is.null(hierDelta)){
      whierDelta <- data.tree::Clone(hierDelta)
    } else whierDelta <- NULL
    
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
          whierBeta[[paste0("level",j)]][[jj]]$betaCons <- betaCons[[paste0("level",j)]][[jj]]$betaCons
          whierBeta[[paste0("level",j)]][[jj]]$workBounds <- workBounds$beta[[paste0("level",j)]][[jj]]$workBounds
          whierBeta[[paste0("level",j)]][[jj]]$fixPar <- fixPar$beta[[paste0("level",j)]][[jj]]$fixPar
          
          whierDelta[[paste0("level",j)]][[jj]]$deltaCons <- deltaCons[[paste0("level",j)]][[jj]]$deltaCons
          whierDelta[[paste0("level",j)]][[jj]]$workBounds <- workBounds$delta[[paste0("level",j)]][[jj]]$workBounds
          whierDelta[[paste0("level",j)]][[jj]]$fixPar <- fixPar$delta[[paste0("level",j)]][[jj]]$fixPar
        }
      } else {
        whierBeta[[paste0("level",j)]]$betaCons <- betaCons[[paste0("level",j)]]$betaCons
        whierBeta[[paste0("level",j)]]$workBounds <- workBounds$beta[[paste0("level",j)]]$workBounds
        whierBeta[[paste0("level",j)]]$fixPar <- fixPar$beta[[paste0("level",j)]]$fixPar
        
        whierDelta[[paste0("level",j)]]$deltaCons <- deltaCons[[paste0("level",j)]]$deltaCons
        whierDelta[[paste0("level",j)]]$workBounds <- workBounds$delta[[paste0("level",j)]]$workBounds
        whierDelta[[paste0("level",j)]]$fixPar <- fixPar$delta[[paste0("level",j)]]$fixPar
      }
    }
    
    cons <- mapCons(whierBeta,whierDelta,hFixPar,hBetaCons,hDeltaCons,hierStates,newformula,formulaDelta,data,mixtures)
    fix  <- mapPar(whierBeta,whierDelta,hFixPar,hBetaCons,cons$deltaCons,hierStates,newformula,formulaDelta,data,mixtures,field="fixPar",check=FALSE)
    par  <- mapPar(whierBeta,whierDelta,hFixPar,hBetaCons,cons$deltaCons,hierStates,newformula,formulaDelta,data,mixtures,field="beta")
    wb   <- mapBounds(whierBeta,whierDelta,hFixPar,hBetaCons,cons$deltaCons,hierStates,newformula,formulaDelta,data,mixtures)
    
    if(!is.list(fixPar)){
      fixPar <- fix
    } else {
      fixPar$beta <- fix$beta
      fixPar$delta <- fix$delta
    }
    
    beta0 <- par$beta
    delta0 <- par$delta
    
    if(!is.list(workBounds)){
      workBounds <- wb
      if(is.null(workBounds$beta)){
        if(any(is.finite(betaLower))){
          workBounds$beta <- cbind(c(betaLower),c(betaUpper))
        }
      }
    } else {
      workBounds$beta <- wb$beta
      workBounds$delta <- wb$delta
    }
    
    if(!is.null(beta0)) {
      if(is.null(whierDelta)) beta0[which(is.na(beta0))] <- 0
    }
    
    if(mixtures==1){
      if(!is.null(hFixPar$beta)) rownames(hFixPar$beta) <- colnames(covs)
      if(!is.null(hFixPar$delta)) rownames(hFixPar$delta) <- colnames(covsDelta)
      if(!is.null(fixPar$beta)) rownames(fixPar$beta) <- colnames(covs)
      if(!is.null(fixPar$delta)) rownames(fixPar$delta) <- colnames(covsDelta)
      if(!is.null(hBetaCons)) rownames(hBetaCons) <- colnames(covs)
      if(!is.null(hDeltaCons)) rownames(hDeltaCons) <- colnames(covsDelta)
      if(!is.null(cons$betaCons)) rownames(cons$betaCons) <- colnames(covs)
      if(!is.null(cons$deltaCons)) rownames(cons$deltaCons) <- colnames(covsDelta)
      if(!is.null(beta0)) rownames(beta0) <- colnames(covs)
      if(!is.null(delta0)) rownames(delta0) <- colnames(covsDelta)
    }
    
    if(mixtures>1 | !is.null(recharge)){
      beta0 <- list(beta=beta0)
      if(mixtures>1) beta0[["pi"]] <- Pi
      if(!is.null(recharge)){
        beta0$g0 <- g0
        beta0$theta <- theta
      }
    }
    
    # populate hierBeta and hierDelta if not provided
    hier <- mapHier(beta0,Pi,delta0,whierBeta,whierDelta,hFixPar,hBetaCons,hDeltaCons,hierStates,newformula,formulaDelta,data,mixtures,recharge)
    whierBeta <- hier$hierBeta
    whierDelta <- hier$hierDelta
    
  } else {
    data <- newForm$data
    whierBeta <- hierBeta
    whierDelta <- hierDelta
    newformula <- hBetaCons <- hDeltaCons <- hFixPar <- workBounds <- cons <- fixPar <- NULL
  }
  if(all(unlist(lapply(recharge,is.null)))) recharge <- NULL

  return(list(data=data,nbStates=nbStates,dist=dist,formula=formula,newformula=newformula,formulaDelta=formulaDelta,beta=beta0,delta=delta0,hierBeta=whierBeta,hierDelta=whierDelta,betaRef=betaRef,betaCons=cons$betaCons,deltaCons=cons$deltaCons,fixPar=fixPar,workBounds=workBounds,stateNames=stateNames,hBetaCons=hBetaCons,hDeltaCons=hDeltaCons,hFixPar=hFixPar,recharge=recharge))

}

nbHierStates <- function(hierStates){
  
  if(!inherits(hierStates,"Node")) stop("'hierStates' must be of class Node; see ?data.tree::Node")
  if(!("state" %in% hierStates$attributesAll)) stop("'hierStates' must include a 'state' field")
  
  hdf <- data.tree::ToDataFrameTypeCol(hierStates, "state")
  if(any(is.na(hdf$state))) stop("'state' field in 'hierStates' cannot contain NAs")
  
  nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))
  
  if(any(duplicated(hdf$state))) stop("'state' field in 'hierStates' cannot contain duplicates")
  if(any(sort(hdf$state)!=1:nbStates)) stop("'state' field in 'hierStates' must include all integers between 1 and ",nbStates)
  if(!data.tree::AreNamesUnique(hierStates)) stop("node names in 'hierStates' must be unique")
  #if(any(is.na(hdf))) stop("missing levels are not permitted in 'hierStates'")
  
  stateNames <- unname(hierStates$Get("name",filterFun=data.tree::isLeaf)[hdf$state])
  
  if(hierStates$height<=2) stop("'hierStates' must contain at least 2 levels below root (i.e., hierStates$height must be > 2)")
  
  if(any(unlist(lapply(data.tree::Traverse(hierStates,traversal="level",filterFun=function(x) !data.tree::isLeaf(x)),function(x) x$count))<2)) stop("each node in 'hierStates' must have at least 2 children")
  
  list(nbStates=nbStates,stateNames=stateNames)
}

getHierDist <- function(hierDist,data,checkData){
  
  if(!inherits(hierDist,"Node")) stop("'hierDist' must be of class Node; see ?data.tree::Node")
  if(!("dist" %in% hierDist$attributesAll)) stop("'hierDist' must include a 'dist' field")
  if(!data.tree::AreNamesUnique(hierDist)) stop("node names in 'hierDist' must be unique")
  if(hierDist$height!=3) stop("'hierDist' hierarchy must contain 2 levels below root (i.e., hierDist$height must be 3)")
  if(!all(hierDist$Get("name",filterFun=function(x) x$level==2)==paste0("level",1:hierDist$count))) stop("hierDist level names from top to bottom should be ",paste0("'level",paste0(1:hierDist$count,"'"),collapse=", ")," (not ",paste0(paste0("'",hierDist$Get("name",filterFun=function(x) x$level==2),"'"),collapse=", "),")")
  nbLevels <- 2*hierDist$count - 1 #hierStates$height #ncol(data.tree::ToDataFrameTypeCol(hierStates))
  
  if(checkData){
    if(is.null(data$level)) stop('data$level must be specified')
    if(!is.factor(data$level)) stop('data$level must be a factor')
    if(nlevels(data$level)!=nbLevels) stop('data$level must contain ',nbLevels,' levels')
    for(k in 1:min(nbLevels,nrow(data))){
      if(data$level[k]!=levels(data$level)[k]) stop("data$level factor levels not ordered correctly; observation ",k," is level ",data$level[k]," but factor level is ",levels(data$level)[k])
    }
    
    ## add additional levels for initial distributions of finer scale states
    #hlevels <- levels(data$level)
    #levels(data$level) <- c(hlevels,paste0("i",hlevels[-1]))
    #for(k in 1:(length(hlevels)-1)){
    #  rInd <- 0
    #  for(j in which(data$level==hlevels[k])){
    #    tmp <- data[j+rInd,]
    #    tmp[names(dist)] <- NA
    #    tmp$level <- factor(paste0("i",hlevels[k+1]),levels=levels(data$level))
    #    data <- DataCombine::InsertRow(data,tmp,j+rInd+1)
    #    rInd <- rInd + 1
    #  }
    #}
    
    if(!all(hierDist$Get("name",filterFun=function(x) x$level==2) %in% paste0("level",levels(data$level)[seq(1,nlevels(data$level),2)]))) 
      stop("'hierDist' level types can only include ",paste(paste0("level",levels(data$level)[seq(1,nlevels(data$level),2)]),collapse=", "))
    
    for(j in gsub("level","",hierDist$Get("name",filterFun=function(x) x$level==2))){
      jDist <- hierDist[[paste0("level",j)]]$Get("dist",filterFun=data.tree::isLeaf)
      if(any(!is.na(jDist))){
        dist <- jDist[!is.na(jDist)]
        distnames <- tmpdistnames <- names(dist)
        if(!all(distnames %in% names(data))){
          for(i in which(is.na(match(distnames,names(data))))){
            if(dist[i] %in% mvndists){
              if(dist[i] %in% c("mvnorm2","rw_mvnorm2")){
                tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"))
              } else if(jDist[i] %in% c("mvnorm3","rw_mvnorm3")){
                tmpdistnames <- c(tmpdistnames[-1],paste0(distnames[i],".x"),paste0(distnames[i],".y"),paste0(distnames[i],".z"))          
              }
            }
          }
          if(any(is.na(match(tmpdistnames,names(data))))) stop(paste0(tmpdistnames[is.na(match(tmpdistnames,names(data)))],collapse=", ")," not found in data")
        }
        for(k in levels(data$level)[-which(levels(data$level)==j)]){
          if(any(!is.na(data[which(data$level==k),tmpdistnames]))) stop(paste(tmpdistnames,collapse=", ")," must be NA for level ",k)
        }
      }
    }
  }# else {
  # lInd <- gsub("level","",sort(hierDist$Get("name",filterFun=function(x) x$level==2),decreasing=TRUE)[1])
  # lLevels <- c(paste0(lInd,"i"),lInd)
  # for(j in sort(hierDist$Get("name",filterFun=function(x) x$level==2),decreasing=TRUE)[-1]){
  #   lLevels <- c(gsub("level","",j),lLevels)
  #   if(j!="level1") {
  #     lLevels <- c(paste0(gsub("level","",j),"i"),lLevels)
  #   }
  # }
  # data <- data.frame(level=factor(lLevels,levels=lLevels))
  #}
  
  dist <- as.list(hierDist$Get("dist",filterFun=data.tree::isLeaf))
  dist <- dist[which(unlist(lapply(dist,function(x) !is.na(x))))]
  dist
}

checkHierFormula <- function(data,hierFormula,hierStates,hierDist,checkData,what="formula"){
  if(is.null(hierFormula)){
    whierFormula <- data.tree::Node$new(hierStates$Get("name",filterFun=data.tree::isRoot))
    whierFormula$AddChild(hierDist$Get("name",filterFun=function(x) x$level==2)[1])
    whierFormula[[hierDist$Get("name",filterFun=function(x) x$level==2)[1]]][[what]] <- ~1
    for(j in hierDist$Get("name",filterFun=function(x) x$level==2)[-1]){
      #if(what=="formula") whierFormula$AddChild(paste0(j,"i"),formula=~1)
      #else {
        whierFormula$AddChild(j)
        whierFormula[[j]][[what]] <- ~1
      #}
    }
  } else {
    if(!inherits(hierFormula,"Node")) stop(ifelse(what=="formula","'hierFormula'","'hierFormulaDelta'")," must be of class Node; see ?data.tree::Node")
    if(!(what %in% hierFormula$attributesAll)) stop(ifelse(what=="formula","'hierFormula'","'hierFormulaDelta'")," must include a ",ifelse(what=="formula","'formula'","'formulaDelta'")," field")
    if(!data.tree::AreNamesUnique(hierFormula)) stop("node names in ",ifelse(what=="formula","'hierFormula'","'hierFormulaDelta'")," must be unique")
    if(hierFormula$height!=2) stop(ifelse(what=="formula","'hierFormula'","'hierFormulaDelta'")," hierarchy must contain 1 level (i.e., ",ifelse(what=="formula","hierFormula","hierFormulaDelta"),"$height must be 2)")
    
    if(checkData){
      if(!all(hierFormula$Get("name",filterFun=function(x) x$level==2) %in% paste0("level",levels(data$level)))) 
        stop("'hierFormula' level types can only include ",paste(paste0("level",levels(data$level)),collapse=", "))
    }
    
    whierFormula <- data.tree::Clone(hierFormula)
    
    #if(!all(hierDist$Get("name",filterFun=function(x) x$level==2)==whierFormula$Get("name",filterFun=function(x) x$level==2)[seq(1,whierFormula$count,2)])) stop("'hierDist' and 'hierFormula' are not consistent; check number of nodes, node names, and node order")
    for(j in hierDist$Get("name",filterFun=function(x) x$level==2)){
      if(is.null(whierFormula[[j]])) {
        whierFormula$AddChild(j)
        whierFormula[[j]][[what]] <- ~1
      } else if(is.null(whierFormula[[j]][[what]])){
        whierFormula[[j]][[what]] <- ~1
      } else if(!is.formula(whierFormula[[j]][[what]])) stop(ifelse(what=="formula","'hierFormula$","'hierFormulaDelta$"),j,ifelse(what=="formula","formula'","formulaDelta'")," must be a formula")
      
      #if(j!=hierDist$Get("name",filterFun=function(x) x$level==2)[1]){
      #  if(is.null(whierFormula[[paste0(j,"i")]])) {
      #    whierFormula$AddChild(paste0(j,"i"),formula=~1)
      #  } else if(is.null(whierFormula[[paste0(j,"i")]]$formula)){
      #    whierFormula[[paste0(j,"i")]]$formula <- ~1
      #  } else if(!is.formula(whierFormula[[paste0(j,"i")]]$formula)) stop("'hierFormula$",paste0(j,"i"),"$formula' must be a formula")
      #}
    }
  }
  whierFormula
}

checkField <- function(what,field,level,hierBeta,hierStates,betaRef,nbCovs,mixtures,initial=FALSE,bounds=FALSE,check=TRUE){
  whierBeta <- data.tree::Clone(hierBeta)
  if(level>1){
    #if(whierBeta[[paste0("level",level)]]$count!=hierStates$Get("count",filterFun = function(x) x$level==(level-1))) stop(what,"$",field," for level",level," must consist of ",hierStates$Get("count",filterFun = function(x) x$level==(level-1))," children: ",paste0(hierStates$Get("name",filterFun=function(x) x$level==level),collapse=", "))
    count <- 0
    for(jj in hierStates$Get("name",filterFun=function(x) x$level==level)){
      inits <- whierBeta[[paste0("level",level)]][[jj]][[field]]
      nlStates <- hierStates$Get("count",filterFun=function(x) x$level==level)[jj]
      jRef <- sum(!(rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==level)[jj],times=hierStates$Get("leafCount",filterFun=function(x) x$level==level)[jj]) %in% betaRef))
      pCount <- ifelse(initial,(nlStates-1),nlStates*(nlStates-1))
      if(is.null(inits) & nlStates) {
        if(check) stop(what,"$level",level,"$",jj,"$",field," are missing")
        else {
          if(field=="fixPar") whierBeta[[paste0("level",level)]][[jj]][[field]] <- inits <- matrix(NA,nbCovs*mixtures,pCount*ifelse(jRef,jRef,1))
          else if(!bounds) whierBeta[[paste0("level",level)]][[jj]][[field]] <- inits <- matrix(count+1:(nbCovs*mixtures*pCount*ifelse(jRef,jRef,1)),nbCovs*mixtures,pCount*ifelse(jRef,jRef,1))
          else whierBeta[[paste0("level",level)]][[jj]][[field]] <- inits <- matrix(c(-Inf,Inf),nbCovs*mixtures*pCount*ifelse(jRef,jRef,1),2,byrow=TRUE)
        }
      } else if(!is.null(inits) & nlStates & field %in% c("betaCons","deltaCons")) {
        if(any(abs(as.integer(inits)-inits)!=0)) stop(what,"$level",level,"$",jj,"$",field," must be a matrix composed of integers")
        if(min(inits)<1 | max(inits)>(nbCovs*mixtures*pCount*ifelse(jRef,jRef,1))) stop(what,"$level",level,"$",jj,"$",field," must be composed of integers between 1 and ",nbCovs*mixtures*pCount*ifelse(jRef,jRef,1))
        whierBeta[[paste0("level",level)]][[jj]][[field]] <- inits <- count + whierBeta[[paste0("level",level)]][[jj]][[field]]
      }
      if(!bounds) {
        if(nlStates && (!is.matrix(inits) || (ncol(inits)!=pCount*ifelse(jRef,jRef,1) | nrow(inits)!=nbCovs*mixtures))) stop(what,"$level",level,"$",jj,"$",field," should consist of ",nbCovs*mixtures," rows and ",pCount*ifelse(jRef,jRef,1), " columns")
      } else {
        if(nlStates && (!is.matrix(inits) || (nrow(inits)!=pCount*ifelse(jRef,jRef,1)*nbCovs*mixtures | ncol(inits)!=2))) stop(what,"$level",level,"$",jj,"$",field," should consist of ",pCount*ifelse(jRef,jRef,1)*nbCovs*mixtures," rows and 2 columns")
      }
      if(nlStates) count <- count + nbCovs*mixtures*pCount*ifelse(jRef,jRef,1)
    }
    inits <- whierBeta[[paste0("level",level)]]$Get(field,filterFun=data.tree::isLeaf,simplify=FALSE)[hierStates$Get("name",filterFun=function(x) x$level==level)]
    if(!bounds) inits <- unlist(inits)
  } else {
    inits <- whierBeta[[paste0("level",level)]][[field]]
    pCount <- ifelse(initial,(hierStates$count-1),hierStates$count*(hierStates$count-1))
    if(is.null(inits)) {
      if(check) stop(what,"$",field," are missing for level",level)
      else {
        if(field=="fixPar") whierBeta[[paste0("level",level)]][[field]] <- inits <- matrix(NA,nbCovs*mixtures,pCount)
        else if(!bounds) whierBeta[[paste0("level",level)]][[field]] <- inits <- matrix(1:(nbCovs*mixtures*pCount),nbCovs*mixtures,pCount)
        else whierBeta[[paste0("level",level)]][[field]] <- inits <- matrix(c(-Inf,Inf),nbCovs*mixtures*pCount,2,byrow=TRUE)
      }
    } else if(field %in% c("betaCons","deltaCons")){
      if(any(abs(as.integer(inits)-inits)!=0)) stop(what,"$level",level,"$",field," must be a matrix composed of integers")
      if(min(inits)<1 | max(inits)>(nbCovs*mixtures*pCount)) stop(what,"$level",level,"$",field," must be composed of integers between 1 and ",nbCovs*mixtures*pCount)
    }
    if(!bounds) {
      if(!is.matrix(inits) || (ncol(inits)!=pCount | nrow(inits)!=nbCovs*mixtures)) stop(what,"$level",level,"$",field," should consist of ",nbCovs*mixtures," rows and ",pCount, " columns")
    } else {
      if(!is.matrix(inits) || (nrow(inits)!=pCount*nbCovs*mixtures | ncol(inits)!=2)) stop(what,"$level",level,"$",field," should consist of ",pCount*nbCovs*mixtures," rows and 2 columns")
    }
    inits <- whierBeta[[paste0("level",level)]][[field]]  
  }
  inits
}

mapCons <- function(hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,check=FALSE){
  
  bc <- betaCons
  dc <- deltaCons
  
  betaRef <- rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  if("betaCons" %in% hierBeta$attributesAll){
    what <- "hierBeta"
    field <- "betaCons"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(stats::model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)

      inits <- checkField(what,field,j,hierBeta,hierStates,betaRef,nbCovs,mixtures,check=check)
      betaInd <- betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]
      naInd <- is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),])
      bc[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][naInd] <- unique(betaInd[naInd])[inits][match(betaInd[naInd],unique(betaInd[naInd]))]

    }
  }
  if("deltaCons" %in% hierDelta$attributesAll){
    what <- "hierDelta"
    field <- "deltaCons"
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        covNames <- colnames(stats::model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
        covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,check=check)
        betaInd <- betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]
        naInd <- is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),])
        bc[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][naInd] <- unique(betaInd[naInd])[inits][match(betaInd[naInd],unique(betaInd[naInd]))]
      } else if(j==1){
        covNames <- colnames(stats::model.matrix(formulaDelta,data))
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,check=check)
        initsInd <- which(is.na(fixPar$delta))#[(mix-1)*nbCovs+1:nbCovs,,drop=FALSE]))
        dc[initsInd] <- initsInd[inits]
      }
    }
  }
  return(list(betaCons=bc, deltaCons = dc))
}

mapPar <- function(hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,field="beta",check=TRUE){
  
  match.arg(field,c("beta","fixPar"))
  
  betaRef <- rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  if(field=="beta") beta <- delta <- NULL
  else {
    beta <- fixPar$beta
    delta <- fixPar$delta
  }
  if(field %in% hierBeta$attributesAll){
    beta <- fixPar$beta
    what <- "hierBeta"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(stats::model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)
      
      inits <- checkField(what,field,j,hierBeta,hierStates,betaRef,nbCovs,mixtures,check=check)
      initsInd <- unique(betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][which(is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]))])
      beta[initsInd] <- inits
    }
  }
  field <- ifelse(field=="beta","delta",field)
  if(field %in% hierDelta$attributesAll){
    if(is.null(beta)) beta <- fixPar$beta
    delta <- fixPar$delta
    what <- "hierDelta"
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        covNames <- colnames(stats::model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
        covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,check=check)
        initsInd <- unique(betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][which(is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]))])
        beta[initsInd] <- inits
      } else if(j==1){
        covNames <- colnames(stats::model.matrix(formulaDelta,data))
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,check=check)
        initsInd <- which(is.na(fixPar$delta))#[(mix-1)*nbCovs+1:nbCovs,,drop=FALSE]))
        delta[initsInd] <- inits
      }
    }
  }
  if(!is.null(beta)) {
    # set any unspecified inital t.p.m. values to zero
    if(field=="delta") beta[is.na(beta)] <- 0
    beta <- matrix(beta[betaCons],nrow(beta),ncol(beta),dimnames=dimnames(betaCons))
  }
  if(!is.null(delta)) delta <- matrix(delta[deltaCons],nrow(delta),ncol(delta),dimnames = dimnames(deltaCons))
  return(list(beta = beta, delta = delta))
}

mapBounds <- function(hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,field="workBounds",check=FALSE){
  
  match.arg(field,"workBounds")
  
  betaRef <- rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  beta <- delta <- NULL
  
  if(field %in% hierBeta$attributesAll){
    betaLower <- matrix(-Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
    betaUpper <- matrix(Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
    what <- "hierBeta"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(stats::model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)
      
      inits <- checkField(what,field,j,hierBeta,hierStates,betaRef,nbCovs,mixtures,bounds=TRUE,check=check)
      initsInd <- unique(betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][which(is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]))])
      if(j>1){
        betaLower[initsInd] <- unlist(lapply(inits,function(x) x[,1]))
        betaUpper[initsInd] <- unlist(lapply(inits,function(x) x[,2]))
      } else {
        betaLower[initsInd] <- inits[,1]
        betaUpper[initsInd] <- inits[,2]
      }
    }
    beta <- cbind(betaLower[betaCons],betaUpper[betaCons])
  }

  if(field %in% hierDelta$attributesAll){
    if(is.null(beta)) {
      betaLower <- matrix(-Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
      betaUpper <- matrix(Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
    }
    deltaLower <- matrix(-Inf,nrow(fixPar$delta),ncol(fixPar$delta),dimnames = dimnames(fixPar$delta))
    deltaUpper <- matrix(Inf,nrow(fixPar$delta),ncol(fixPar$delta),dimnames = dimnames(fixPar$delta))
    
    what <- "hierDelta"
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        covNames <- colnames(stats::model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
        covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
        nbCovs <- length(covNames)

        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,bounds=TRUE,check=check)
        initsInd <- unique(betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][which(is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]))])
        betaLower[initsInd] <- unlist(lapply(inits,function(x) x[,1]))
        betaUpper[initsInd] <- unlist(lapply(inits,function(x) x[,2]))
      } else if(j==1){
        covNames <- colnames(stats::model.matrix(formulaDelta,data))
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,bounds=TRUE,check=check)
        initsInd <- which(is.na(fixPar$delta))#[(mix-1)*nbCovs+1:nbCovs,,drop=FALSE]))
        deltaLower[initsInd] <- inits[,1]
        deltaUpper[initsInd] <- inits[,2]
      }
    }
    beta <- cbind(betaLower[betaCons],betaUpper[betaCons])
    delta <- cbind(deltaLower[deltaCons],deltaUpper[deltaCons])
  }
  return(list(beta = beta, delta = delta))
}

getRechargeParms <- function(recharge,data,hierBeta){

  hierG0 <- hierTheta <- list()
  g0 <- numeric()
  theta <- numeric()
  for(j in names(recharge)){
    tmpg0covs <- stats::model.matrix(recharge[[j]]$g0,data)
    tmpnbG0 <- ncol(tmpg0covs)
    tmpG0 <- hierBeta[[j]]$g0
    if(!is.null(tmpG0)){
      if(length(tmpG0)!=tmpnbG0 | any(!is.numeric(tmpG0))) stop("hierBeta$",j,"$g0 must be a numeric vector of length ",tmpnbG0)
      hierG0[[j]] <- tmpG0
    } else {
      hierG0[[j]] <- rep(0,tmpnbG0)
    }
    names(hierG0[[j]]) <- colnames(tmpg0covs)
    g0 <- c(g0,hierG0[[j]])
    
    tmprecovs <- stats::model.matrix(recharge[[j]]$theta,data)
    tmpnbRecovs <- ncol(tmprecovs)
    tmpTheta <- hierBeta[[j]]$theta
    if(!is.null(tmpTheta)){
      if(length(tmpTheta)!=tmpnbRecovs | any(!is.numeric(tmpTheta))) stop("hierBeta$",j,"$theta must be a numeric vector of length ",tmpnbRecovs)
      hierTheta[[j]] <- tmpTheta
    } else {
      hierTheta[[j]] <- c(-1,rep(0,tmpnbRecovs-1))
    }
    names(hierTheta[[j]]) <- colnames(tmprecovs)
    theta <- c(theta,hierTheta[[j]])
  }

  list(g0=g0,theta=theta,hierG0=hierG0,hierTheta=hierTheta)
}

mapHier <- function(beta,Pi,delta,hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,recharge=NULL,fill=FALSE){
  
  if(is.null(hierBeta)){
    hierBeta <- data.tree::Node$new("hierBeta")
    hierBeta$AddChild(paste0("level1"))
    for(j in 2:(hierStates$height-1)){
      hierBeta$AddChild(paste0("level",j))
      for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
        hierBeta[[paste0("level",j)]]$AddChild(jj)
      }
    }
  }
  if(is.null(hierDelta)){
    hierDelta <- data.tree::Node$new("hierDelta")
    hierDelta$AddChild(paste0("level1"))
    for(j in 2:(hierStates$height-1)){
      hierDelta$AddChild(paste0("level",j))
      for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
        hierDelta[[paste0("level",j)]]$AddChild(jj)
      }
    }        
  }
  
  if(is.null(beta)){
    beta <- matrix(0,nrow(fixPar$beta),ncol(fixPar$beta),dimnames=dimnames(fixPar$beta))
  } else if(is.list(beta)){
    if(is.null(beta$beta)) beta <- matrix(0,nrow(fixPar$beta),ncol(fixPar$beta),dimnames=dimnames(fixPar$beta))
    else beta <- beta$beta
  }  
  if(is.null(delta)) delta <- matrix(0,nrow(fixPar$delta),ncol(fixPar$delta),dimnames=dimnames(fixPar$delta))
  
  betaRef <- rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  whierBeta <- NULL
  if(is.list(hierBeta)){
    if(inherits(hierBeta$beta,"Node")) whierBeta <- data.tree::Clone(hierBeta$beta)
  } else if(inherits(hierBeta,"Node")) whierBeta <- data.tree::Clone(hierBeta)

  if(inherits(whierBeta,"Node")){
    what <- "hierBeta"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(stats::model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)
      if(mixtures>1) covNames <- paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs))
      if(j>1){
        initsInd <- betaCons[covNames,][which(is.na(fixPar$beta[covNames,]))]
        dimNames <- mapply(function(x) mapply(`[[`, dimnames(fixPar$beta), arrayInd(x, dim(fixPar$beta))),initsInd)
        inits <- beta[initsInd]
        count <- 0
        t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
        names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
        for(jj in names(t)){
          nStates <- t[[jj]]$count#length(t[[jj]]$Get("state",filterFun = data.tree::isLeaf))
          jRef <- sum(!(rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==j)[jj],times=hierStates$Get("leafCount",filterFun=function(x) x$level==j)[jj]) %in% betaRef))
          pCount <- nStates*(nStates-1)*ifelse(jRef,jRef,1)
          if(nStates){
            if(is.null(whierBeta[[paste0("level",j)]][[jj]]$betaCons)) {
              iRef <- count + matrix(1:(pCount*nbCovs*mixtures),nbCovs*mixtures,pCount)
            } else iRef <- count + whierBeta[[paste0("level",j)]][[jj]]$betaCons
            idimNames <- list(covNames,unique(dimNames[2,count + matrix(1:(pCount*nbCovs*mixtures),nbCovs*mixtures,pCount)]))
            if(fill & is.null(whierBeta[[paste0("level",j)]][[jj]]$betaCons)) whierBeta[[paste0("level",j)]][[jj]]$betaCons <- matrix(iRef-count,nbCovs*mixtures,pCount,dimnames=idimNames) 
            if(fill & is.null(whierBeta[[paste0("level",j)]][[jj]]$fixPar)) whierBeta[[paste0("level",j)]][[jj]]$fixPar <- matrix(NA,nbCovs*mixtures,pCount,dimnames=idimNames) 
            if(fill & is.null(whierBeta[[paste0("level",j)]][[jj]]$workBounds)) whierBeta[[paste0("level",j)]][[jj]]$workBounds <- matrix(c(-Inf,Inf),nbCovs*mixtures*pCount,2,byrow=TRUE, dimnames = list(paste0(idimNames[[1]],":",rep(idimNames[[2]],each=length(idimNames[[1]]))),c("lower","upper")))
            betaInits <- inits[iRef]
            if(!is.null(whierBeta[[paste0("level",j)]][[jj]]$betaCons)) betaInits <- betaInits[whierBeta[[paste0("level",j)]][[jj]]$betaCons]
            whierBeta[[paste0("level",j)]][[jj]]$beta <- matrix(betaInits,nbCovs*mixtures,pCount,dimnames=idimNames)
            count <- count + (pCount*nbCovs*mixtures)
          } else {
            if(is.null(whierBeta[[paste0("level",j)]][[jj]])){
              #if(fill) whierBeta[[paste0("level",j)]]$AddChild(jj,betaCons=NULL)
            } else if(!is.null(whierBeta[[paste0("level",j)]][[jj]]$betaCons)) stop("There should be no parameters for hierBeta$level",j,"$",jj,"$betaCons")
            if(is.null(whierBeta[[paste0("level",j)]][[jj]])){
              #if(fill) whierBeta[[paste0("level",j)]]$AddChild(jj,fixPar=NULL)
            } else if(!is.null(whierBeta[[paste0("level",j)]][[jj]]$fixPar)) stop("There should be no parameters for hierBeta$level",j,"$",jj,"$fixPar")
            if(is.null(whierBeta[[paste0("level",j)]][[jj]])){
              #if(fill) whierBeta[[paste0("level",j)]]$AddChild(jj,workBounds=NULL)
            } else if(!is.null(whierBeta[[paste0("level",j)]][[jj]]$workBounds)) stop("There should be no parameters for hierBeta$level",j,"$",jj,"$workBounds")
            if(is.null(whierBeta[[paste0("level",j)]][[jj]])){
              #whierBeta[[paste0("level",j)]]$AddChild(jj,beta=NULL)
            } else if(!is.null(whierBeta[[paste0("level",j)]][[jj]]$beta)) stop("There should be no parameters for hierBeta$level",j,"$",jj,"$beta")
          }
        }
      } else {
        initsInd <- unique(betaCons[covNames,][which(is.na(fixPar$beta[covNames,]))])
        dimNames <- mapply(function(x) mapply(`[[`, dimnames(fixPar$beta), arrayInd(x, dim(fixPar$beta))),initsInd)
        dimNames <- list(covNames,unique(dimNames[2,]))
        betaInits <- beta[initsInd]
        if(!is.null(whierBeta[[paste0("level",j)]]$betaCons)) betaInits <- betaInits[whierBeta[[paste0("level",j)]]$betaCons]
        if(fill & is.null(whierBeta[[paste0("level",j)]]$betaCons)) whierBeta[[paste0("level",j)]]$betaCons <- matrix(1:length(betaInits),nbCovs*mixtures,dimnames = dimNames)
        if(fill & is.null(whierBeta[[paste0("level",j)]]$fixPar)) whierBeta[[paste0("level",j)]]$fixPar <- matrix(rep(NA,length(betaInits)),nbCovs*mixtures,dimnames = dimNames)
        if(fill & is.null(whierBeta[[paste0("level",j)]]$workBounds)) whierBeta[[paste0("level",j)]]$workBounds <- matrix(c(-Inf,Inf),length(betaInits),2,byrow=TRUE, dimnames = list(paste0(dimNames[[1]],":",rep(dimNames[[2]],each=length(dimNames[[1]]))),c("lower","upper")))
        whierBeta[[paste0("level",j)]]$beta <- matrix(betaInits,nbCovs*mixtures,dimnames = dimNames)
      }
    }
  }
  
  whierDelta <- NULL
  if(inherits(hierDelta,"Node")){
    whierDelta <- data.tree::Clone(hierDelta)
    what <- "hierDelta"
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        covNames <- colnames(stats::model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
        covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
        nbCovs <- length(covNames)
        if(mixtures>1) covNames <- paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs))
        initsInd <- betaCons[covNames,][which(is.na(fixPar$beta[covNames,]))]
        dimNames <- mapply(function(x) mapply(`[[`, dimnames(fixPar$beta), arrayInd(x, dim(fixPar$beta))),initsInd)
        inits <- beta[initsInd]
        count <- 0
        t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
        names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
        for(jj in names(t)){
          nStates <- t[[jj]]$count#length(t[[jj]]$Get("state",filterFun = data.tree::isLeaf))
          jRef <- sum(!(rep(hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==j)[jj],times=hierStates$Get("leafCount",filterFun=function(x) x$level==j)[jj]) %in% betaRef))
          pCount <- (nStates-1)*ifelse(jRef,jRef,1)
          if(nStates){
            if(is.null(whierDelta[[paste0("level",j)]][[jj]]$deltaCons)) {
              iRef <- count + matrix(1:(pCount*nbCovs*mixtures),nbCovs*mixtures,pCount)
            } else iRef <- count + whierDelta[[paste0("level",j)]][[jj]]$deltaCons
            idimNames <- list(covNames,gsub('.*->',"state",unique(dimNames[2,count + matrix(1:(pCount*nbCovs*mixtures),nbCovs*mixtures,pCount)])))
            if(fill & is.null(whierDelta[[paste0("level",j)]][[jj]]$deltaCons)) whierDelta[[paste0("level",j)]][[jj]]$deltaCons <- matrix(iRef-count,nbCovs*mixtures,pCount,dimnames=idimNames)
            if(fill & is.null(whierDelta[[paste0("level",j)]][[jj]]$fixPar)) whierDelta[[paste0("level",j)]][[jj]]$fixPar <- matrix(NA,nbCovs*mixtures,pCount,dimnames=idimNames)
            if(fill & is.null(whierDelta[[paste0("level",j)]][[jj]]$workBounds)) whierDelta[[paste0("level",j)]][[jj]]$workBounds <- matrix(c(-Inf,Inf),nbCovs*mixtures*pCount,2,byrow=TRUE, dimnames = list(paste0(idimNames[[1]],":",rep(idimNames[[2]],each=length(idimNames[[1]]))),c("lower","upper")))
            deltaInits <- inits[iRef]
            if(!is.null(whierDelta[[paste0("level",j)]][[jj]]$deltaCons)) deltaInits <- deltaInits[whierDelta[[paste0("level",j)]][[jj]]$deltaCons]
            whierDelta[[paste0("level",j)]][[jj]]$delta <- matrix(deltaInits,nbCovs*mixtures,pCount,dimnames=idimNames)
            count <- count + pCount*nbCovs*mixtures
          } else {
            if(is.null(whierDelta[[paste0("level",j)]][[jj]])){
              #if(fill) whierDelta[[paste0("level",j)]]$AddChild(jj,deltaCons=NULL)
            } else if(!is.null(whierDelta[[paste0("level",j)]][[jj]]$deltaCons)) stop("There should be no parameters for hierDelta$level",j,"$",jj,"$deltaCons")
            if(is.null(whierDelta[[paste0("level",j)]][[jj]])){
              #if(fill) whierDelta[[paste0("level",j)]]$AddChild(jj,fixPar=NULL)
            } else if(!is.null(whierDelta[[paste0("level",j)]][[jj]]$fixPar)) stop("There should be no parameters for hierDelta$level",j,"$",jj,"$fixPar")
            if(is.null(whierDelta[[paste0("level",j)]][[jj]])){
              #if(fill) whierDelta[[paste0("level",j)]]$AddChild(jj,workBounds=NULL)
            } else if(!is.null(whierDelta[[paste0("level",j)]][[jj]]$workBounds)) stop("There should be no parameters for hierDelta$level",j,"$",jj,"$workBounds")
            if(is.null(whierDelta[[paste0("level",j)]][[jj]])){
              #whierDelta[[paste0("level",j)]]$AddChild(jj,delta=NULL)
            } else if(!is.null(whierDelta[[paste0("level",j)]][[jj]]$delta)) stop("There should be no parameters for hierDelta$level",j,"$",jj,"$delta")
          }
        }
      } else if(j==1){
        covNames <- colnames(stats::model.matrix(formulaDelta,data))
        nbCovs <- length(covNames)
        if(mixtures>1) covNames <- paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs))
        initsInd <- unique(deltaCons[covNames,][which(is.na(fixPar$delta[covNames,]))])
        dimNames <- mapply(function(x) mapply(`[[`, dimnames(fixPar$delta), arrayInd(x, dim(fixPar$delta))),initsInd)
        dimNames <- list(covNames,unique(dimNames[2,]))
        deltaInits <- delta[initsInd]
        if(!is.null(whierDelta[[paste0("level",j)]]$deltaCons)) deltaInits <- deltaInits[whierDelta[[paste0("level",j)]]$deltaCons]
        if(fill & is.null(whierDelta[[paste0("level",j)]]$deltaCons)) whierDelta[[paste0("level",j)]]$deltaCons <- matrix(1:length(deltaInits),nbCovs*mixtures,dimnames=dimNames)
        if(fill & is.null(whierDelta[[paste0("level",j)]]$fixPar)) whierDelta[[paste0("level",j)]]$fixPar <- matrix(rep(NA,length(deltaInits)),nbCovs*mixtures,dimnames=dimNames)
        if(fill & is.null(whierDelta[[paste0("level",j)]]$workBounds)) whierDelta[[paste0("level",j)]]$workBounds <- matrix(c(-Inf,Inf),length(deltaInits),2,byrow=TRUE, dimnames = list(paste0(dimNames[[1]],":",rep(dimNames[[2]],each=length(dimNames[[1]]))),c("lower","upper")))
        whierDelta[[paste0("level",j)]]$delta <- matrix(deltaInits,nbCovs*mixtures,dimnames=dimNames)
      }
    }
  }
  
  if(!is.null(recharge)){
    reParms <- getRechargeParms(recharge,data,hierBeta)
    for(j in names(recharge)){
      whierBeta[[j]]$g0 <- reParms$hierG0[[j]]
      whierBeta[[j]]$theta <- reParms$hierTheta[[j]]
    }
  }
  
  if(mixtures>1){
    whierBeta <- list(beta=whierBeta, pi=Pi)
  }
  return(list(hierBeta=whierBeta,hierDelta=whierDelta))
}

find_intercept <- function(intCov,otherCovs,factorCovs,nfInd,data){
  if(!length(intCov)){
    fMat <- matrix(unlist(lapply(factorCovs,function(x) grepl(x,otherCovs))),nrow=length(factorCovs),ncol=length(otherCovs),byrow=TRUE,dimnames=list(factorCovs,otherCovs))
    fCovs <- which(apply(fMat,2,function(x) any(x)))
    if(!length(fCovs)) stop("hierFormula and hierFormulaDelta formulas must contain an intercept (unless the formula contains a factor covariate)")
    nfCovs <- which(apply(matrix(unlist(lapply(nfInd,function(x) grepl(x,otherCovs))),nrow=length(nfInd),ncol=length(otherCovs),byrow=TRUE),2,function(x) any(x)))
    fMat[,otherCovs[nfCovs]] <- FALSE
    if(sum(rowSums(fMat)>0)>1){
      factorTerms <- factorCovs[which(rowSums(fMat)>0)]
      multiFact <- which(colSums(fMat)>1)
      for(k in which(rowSums(fMat)>0)){
        if(sum(fMat[k,])!=length(fCovs)){
          if(sum(fMat[k,colSums(fMat)==1]) && sum(fMat[k,colSums(fMat)==1])!=nlevels(data[[factorCovs[k]]])) {
            fMat[k,] <- FALSE
          } else {
            if(sum(fMat)>sum(fMat[,multiFact])) fMat[k,multiFact] <- FALSE
          }
        }
      }
    }
    fCovs <- which(apply(fMat,2,function(x) any(x)))
    intCov <- otherCovs[fCovs]
    otherCovs <- otherCovs[-fCovs]
  }
  list(intCov=intCov,otherCovs=otherCovs)
}

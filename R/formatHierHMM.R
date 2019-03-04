
#' Convert hierarchical HMM structure to a conventional HMM
#' 
#' @param data \code{\link{momentuHierHMMData}} object or a data frame containing the data streams and covariates. 
#' @param hierStates A hierarchical data structure \code{\link[data.tree]{Node}} for the states ('state').  See \code{\link{fitHMM}}. 
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). See \code{\link{fitHMM}}. 
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of regression coefficients for the transition probabilities at each level of the hierarchy, including initial values ('beta'), parameter equality constraints ('betaCons'), fixed parameters ('fixPar'), and working scale bounds ('workBounds'). See \code{\link{fitHMM}}. 
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the initial distribution at each level of the hierarchy, including initial values ('delta'), parameter equality constraints ('deltaCons'), fixed parameters ('fixPar'), and working scale bounds ('workBounds'). See \code{\link{fitHMM}}. 
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). See \code{\link{fitHMM}}. Default: \code{NULL} (only hierarchical-level effects, with no covariate effects). 
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). See \code{\link{fitHMM}}. Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). See \code{\link{fitHMM}}. Default: \code{mixtures=1}.  
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
#' @importFrom data.tree Node Get Do ToDataFrameTypeCol Traverse Aggregate AreNamesUnique isRoot isLeaf Clone
#' @importFrom stats terms
formatHierHMM <- function(data,hierStates,hierDist,
                          hierBeta=NULL,hierDelta=NULL,
                          hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,
                          checkData=TRUE){
  
  if(is.null(data)) checkData <- FALSE
  
  if(!inherits(hierStates,"Node")) stop("'hierStates' must be of class Node; see ?data.tree::Node")
  if(!("state" %in% hierStates$fieldsAll)) stop("'hierStates' must include a 'state' field")
  
  hdf <- data.tree::ToDataFrameTypeCol(hierStates, "state")
  if(any(is.na(hdf$state))) stop("'state' field in 'hierStates' cannot contain NAs")
  nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))
  
  if(any(duplicated(hdf$state))) stop("'state' field in 'hierStates' cannot contain duplicates")
  if(any(sort(hdf$state)!=1:nbStates)) stop("'state' field in 'hierStates' must include all integers between 1 and ",nbStates)
  if(!data.tree::AreNamesUnique(hierStates)) stop("node names in 'hierStates' must be unique")
  #if(any(is.na(hdf))) stop("missing levels are not permitted in 'hierStates'")
  
  stateNames <- unname(hierStates$Get("name",filterFun=data.tree::isLeaf)[hdf$state])
  
  if(hierStates$height<=2) stop("'hierStates' must contain at least 2 levels below root (i.e., hierStates$height must be > 2)")
  
  if(any(unlist(lapply(Traverse(hierStates,traversal="level",filterFun=function(x) !isLeaf(x)),function(x) x$count))<2)) stop("each node in 'hierStates' must have at least 2 children")
  
  if(!inherits(hierDist,"Node")) stop("'hierDist' must be of class Node; see ?data.tree::Node")
  if(!("dist" %in% hierDist$fieldsAll)) stop("'hierDist' must include a 'dist' field")
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
      for(k in levels(data$level)[-which(levels(data$level)==j)]){
        if(any(!is.na(data[which(data$level==k),names(hierDist[[paste0("level",j)]]$Get("dist",filterFun=isLeaf))]))) stop(paste(names(hierDist[[paste0("level",j)]]$Get("dist",filterFun=isLeaf)),collapse=", ")," must be NA for level ",k)
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

  whierFormulaDelta <- checkHierFormula(data,hierFormulaDelta,hierStates,hierDist,checkData,what="formulaDelta")
  formulaDelta <- whierFormulaDelta$level1$formulaDelta
  whierFormulaDelta$RemoveChild("level1")
  
  whierFormula <- checkHierFormula(data,hierFormula,hierStates,hierDist,checkData)
  for(j in names(whierFormulaDelta$children)){
    whierFormula$AddChild(paste0(j,"i"),formula=whierFormulaDelta[[j]]$formulaDelta)
  }
  formula <- formatHierFormula(data,whierFormula)
  
  recharge <- newFormulas(formula,nbStates)$recharge
  if(!is.null(recharge)) stop("sorry, hierarchical recharge models are not currently supported")
  
  # set t.p.m. reference states based on top level
  betaRef <- rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))

  beta0 <- delta0 <- NULL
  
  if(!is.null(data)){
    
    covs <- model.matrix(formula,data)
    nbCovs <- ncol(covs)
    
    fixPar <- list()
    fixPar$beta <- matrix(NA,nbCovs*mixtures,nbStates*(nbStates-1))
    rownames(fixPar$beta) <- paste0(rep(colnames(covs),mixtures),"_mix",rep(1:mixtures,each=length(colnames(covs))))
    colnames(fixPar$beta) <- c(sapply(1:nbStates,function(x) paste(rep(x,each=nbStates-1),"->",1:nbStates)[-betaRef[x]]))
    
    betaCons <- matrix(1:(nbCovs*mixtures*nbStates*(nbStates-1)),nbCovs*mixtures,nbStates*(nbStates-1))
    dimnames(betaCons)<-dimnames(fixPar$beta)
    
    aInd <- NULL
    nbAnimals <- length(unique(data$ID))
    for(i in 1:nbAnimals)
      aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
    
    covsDelta <- model.matrix(formulaDelta,data[aInd,,drop=FALSE])
    nbCovsDelta <- ncol(covsDelta)
    fixPar$delta <- matrix(NA,nbCovsDelta*mixtures,nbStates-1,byrow=TRUE)
    colnames(fixPar$delta) <- paste("state",2:nbStates)
    rownames(fixPar$delta) <- paste0(rep(colnames(covsDelta),mixtures),"_mix",rep(1:mixtures,each=length(colnames(covsDelta))))
    #if(any(betaRef==1)){
    for(mix in 1:mixtures){
      fixPar$delta[paste0("(Intercept)","_mix",mix),(2:nbStates)[-(betaRef-1)]-1] <- -1.e+10
      if(nbCovsDelta>1) fixPar$delta[paste0(colnames(covsDelta)[which(colnames(covsDelta)=="(Intercept)")],"_mix",mix),(2:nbStates)[-(betaRef-1)]-1] <- 0
    }
    #} else {
    #  if(is.null(workBounds$delta)){
    #    deltaLower <- matrix(-Inf,nbCovsDelta*mixtures,nbStates-1,dimnames=list(rownames(fixPar$delta)))
    #    deltaUpper <- matrix( Inf,nbCovsDelta*mixtures,nbStates-1,dimnames=list(rownames(fixPar$delta)))
    #    for(mix in 1:mixtures){
    #     deltaLower[paste0("level",levels(data$level)[1],"_mix",mix),betaRef-1] <- 100
    #     fixPar$delta[paste0("level",levels(data$level)[1],"_mix",mix),-(betaRef-1)] <- 0
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
    ointCov <- c(paste0("level",levels(data$level)[(2*j):nbLevels]),paste0("I((level == \"",levels(data$level)[(2*j):nbLevels],"\") * 1)"))
    ointCov <- ointCov[which(ointCov %in% colnames(covs))]
    levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-1],"\")"),colnames(covs),fixed=TRUE)]
    levelCovs <- levelCovs[which(!(levelCovs %in% intCov))]
    otherCovs <- colnames(covs)[-which(colnames(covs) %in% c(intCov,levelCovs,ointCov))]
    for(k in names(t)){
      tt <- data.tree::Traverse(t[[k]],filterFun=function(x) x$level==j+1)
      withinConstr <- acrossConstr <- list()
      if(length(tt)){
        for(h in 1:length(tt)){
          levelStates <- tt[[h]]$Get("state",filterFun = data.tree::isLeaf)
          stateInd <- tt[[h]]$Get("levelStates",filterFun=data.tree::isLeaf)
          withinConstr[[h]] <- match(paste0(rep(levelStates,each=length(stateInd))," -> ",stateInd),colnames(fixPar$beta),nomatch=0)
          if(any(withinConstr[[h]])){
            for(mix in 1:mixtures){
              betaCons[paste0(intCov,"_mix",mix),withinConstr[[h]]] <- min(betaCons[paste0(intCov,"_mix",mix),withinConstr[[h]]])
              for(lcovs in levelCovs){
                betaCons[paste0(lcovs,"_mix",mix),withinConstr[[h]]] <- min(betaCons[paste0(lcovs,"_mix",mix),withinConstr[[h]]])
              }
              fixPar$beta[paste0(intCov,"_mix",mix),withinConstr[[h]]] <- -1.e+10
              for(lcovs in levelCovs){
                fixPar$beta[paste0(lcovs,"_mix",mix),withinConstr[[h]]] <- 0
              }
            }
          }
          # constrain across transitions to reference states
          levelStates <- lapply(tt[(1:length(tt))[-h]],function(x) x$Get("state",filterFun = data.tree::isLeaf))
          for(ll in 1:length(levelStates)){
            acrossConstr[[h]] <- match(paste0(rep(levelStates[[ll]],each=length(stateInd[-1]))," -> ",stateInd[-1]),colnames(fixPar$beta),nomatch=0)
            acrossRef <- match(paste0(levelStates[[ll]]," -> ",stateInd[1]),colnames(fixPar$beta),nomatch=0)
            for(mix in 1:mixtures){
              if(any(acrossConstr[[h]])){
                betaCons[paste0(intCov,"_mix",mix),acrossConstr[[h]]] <- min(betaCons[paste0(intCov,"_mix",mix),acrossConstr[[h]]])
                for(lcovs in levelCovs){
                  betaCons[paste0(lcovs,"_mix",mix),acrossConstr[[h]]] <- min(betaCons[paste0(lcovs,"_mix",mix),acrossConstr[[h]]])
                }
                betaCons[paste0(ointCov,"_mix",mix),acrossConstr[[h]]] <- min(betaCons[paste0(ointCov,"_mix",mix),acrossConstr[[h]]])
                for(ocovs in otherCovs){
                  betaCons[paste0(ocovs,"_mix",mix),acrossConstr[[h]]] <- min(betaCons[paste0(ocovs,"_mix",mix),acrossConstr[[h]]])
                }
              }
              if(any(acrossRef)){
                betaCons[paste0(intCov,"_mix",mix),acrossRef] <- min(betaCons[paste0(intCov,"_mix",mix),acrossRef])
                for(lcovs in levelCovs){
                  betaCons[paste0(lcovs,"_mix",mix),acrossRef] <- min(betaCons[paste0(lcovs,"_mix",mix),acrossRef])
                }
                betaCons[paste0(ointCov,"_mix",mix),acrossRef] <- min(betaCons[paste0(ointCov,"_mix",mix),acrossRef])
                for(ocovs in otherCovs){
                  betaCons[paste0(ocovs,"_mix",mix),acrossRef] <- min(betaCons[paste0(ocovs,"_mix",mix),acrossRef])
                }
              }
              if(any(acrossConstr[[h]])){
                fixPar$beta[paste0(intCov,"_mix",mix),acrossConstr[[h]]] <- -1.e+10
                for(lcovs in levelCovs){
                  fixPar$beta[paste0(lcovs,"_mix",mix),acrossConstr[[h]]] <- 0
                }
                fixPar$beta[paste0(ointCov,"_mix",mix),acrossConstr[[h]]] <- -1.e+10
                for(ocovs in otherCovs){
                  fixPar$beta[paste0(ocovs,"_mix",mix),acrossConstr[[h]]] <- 0
                }
              }
              if(any(acrossRef)){
                fixPar$beta[paste0(ointCov,"_mix",mix),acrossRef] <- -1.e+10
                for(ocovs in otherCovs){
                  fixPar$beta[paste0(ocovs,"_mix",mix),acrossRef] <- 0
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
    
    betaLower <- matrix(-Inf,nbCovs*mixtures,nbStates*(nbStates-1),dimnames = dimnames(betaCons))
    betaUpper <- matrix( Inf,nbCovs*mixtures,nbStates*(nbStates-1),dimnames = dimnames(betaCons))
    
    for(j in 2:(hierStates$height-1)){
      
      t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
      names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
      
      intCov <- c(paste0("level",levels(data$level)[2*j-2]),paste0("I((level == \"",levels(data$level)[2*j-2],"\") * 1)"))
      intCov <- intCov[which(intCov %in% colnames(covs))]
      levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-2],"\")"),colnames(covs),fixed=TRUE)]
      levelCovs <- levelCovs[which(!(levelCovs %in% intCov))]
      
      #initial distribution       
      for(k in names(t)){
        levelStates <- t[[k]]$Get("state",filterFun = data.tree::isLeaf)
        allStates <- levelStateFun(t[[k]])
        fromState <- data.tree::Aggregate(t[[k]],"state",min)
        tt <- data.tree::Traverse(t[[k]],filterFun=function(x) x$level==j+1)
        toStates <- unlist(lapply(tt,function(x) data.tree::Aggregate(x,"state",min)))#data.tree::Aggregate(tt[[h]],"state",min)
        if(is.null(toStates)) toStates <- t[[k]]$state
        allTrans <- paste0(rep(levelStates,each=length(allStates[which(!(allStates %in% betaRef[toStates]))]))," -> ",allStates[which(!(allStates %in% betaRef[toStates]))])
        initConstr <- match(allTrans[which(!(allTrans %in% paste0(rep(fromState,each=length(toStates))," -> ",toStates)))],colnames(fixPar$beta),nomatch=0)
        refConstr <- match(allTrans[which(allTrans %in% paste0(rep(fromState,each=length(toStates))," -> ",toStates))],colnames(fixPar$beta),nomatch=0)
        for(mix in 1:mixtures){
          if(any(initConstr)){
            betaCons[paste0(intCov,"_mix",mix),initConstr] <- min(betaCons[paste0(intCov,"_mix",mix),initConstr])
            for(lcovs in levelCovs){
              betaCons[paste0(lcovs,"_mix",mix),initConstr] <- min(betaCons[paste0(lcovs,"_mix",mix),initConstr])
            }
            fixPar$beta[paste0(intCov,"_mix",mix),initConstr] <- -1.e+10
            for(lcovs in levelCovs){
              fixPar$beta[paste0(lcovs,"_mix",mix),initConstr] <- 0
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
              } else if(!length(levelCovs)) fixPar$beta[paste0(intCov,"_mix",mix),refConstr] <- 500
            }
          }
        }
      }
      
      intCov <- c(paste0("level",levels(data$level)[2*j-1]),paste0("I((level == \"",levels(data$level)[2*j-1],"\") * 1)"))
      intCov <- intCov[which(intCov %in% colnames(covs))]
      levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-1],"\")"),colnames(covs),fixed=TRUE)]
      levelCovs <- levelCovs[which(!(levelCovs %in% intCov))]

      # t.p.m.
      for(k in names(t)){  
        levelStates <- t[[k]]$Get("state",filterFun = data.tree::isLeaf)
        allStates <- levelStateFun(t[[k]])
        tt <- data.tree::Traverse(t[[k]],filterFun=function(x) x$level==j+1)
        fromStates <- toStates <- unlist(lapply(tt,function(x) data.tree::Aggregate(x,"state",min)))
        if(is.null(toStates)) toStates <- fromStates <- t[[k]]$state
        allTrans <- paste0(rep(levelStates,each=length(allStates[which(!(allStates %in% betaRef[toStates]))]))," -> ",allStates[which(!(allStates %in% betaRef[toStates]))])
        initConstr <- match(allTrans[which(!(allTrans %in% paste0(rep(fromStates,each=length(toStates))," -> ",toStates)))],colnames(fixPar$beta),nomatch=0)
        refConstr <- match(allTrans[which(allTrans %in% paste0(rep(fromStates,each=length(toStates))," -> ",toStates))],colnames(fixPar$beta),nomatch=0)
        for(mix in 1:mixtures){
          if(any(initConstr)){
            betaCons[paste0(intCov,"_mix",mix),initConstr] <- min(betaCons[paste0(intCov,"_mix",mix),initConstr])
            for(lcovs in levelCovs){
              betaCons[paste0(lcovs,"_mix",mix),initConstr] <- min(betaCons[paste0(lcovs,"_mix",mix),initConstr])
            }
            fixPar$beta[paste0(intCov,"_mix",mix),initConstr] <- -1.e+10
            for(lcovs in levelCovs){
              fixPar$beta[paste0(lcovs,"_mix",mix),initConstr] <- 0
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
              } else if(!length(levelCovs)) fixPar$beta[paste0(intCov,"_mix",mix),refConstr] <- 500
            }
          }
        }
      }
    }
    fixInd <- which(fixPar$beta==-1.e+10)
    if(length(fixInd)) betaCons[fixInd] <- fixInd[1]
    
    Pi <- g0 <- theta <- NULL
    if(!is.null(hierBeta)){
      if(mixtures>1){
        if(!is.list(hierBeta) || !all(names(hierBeta) %in% c("beta","pi","g0","theta"))) stop("hierBeta must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
        Pi <- hierBeta$pi
      }
      if(!is.null(recharge)){
        if(!is.list(hierBeta) || !all(names(hierBeta) %in% c("beta","pi","g0","theta"))) stop("hierBeta must be a list with elements named 'beta', 'g0', and/or 'theta' when including a recharge model")
        g0 <- hierBeta$g0
        theta <- hierBeta$theta
      }
      if(is.list(hierBeta)) hierBeta <- hierBeta$beta
    }
    
    cons <- mapCons(hierBeta,hierDelta,fixPar,betaCons,hierStates,formula,formulaDelta,data,mixtures)
    fix  <- mapPar(hierBeta,hierDelta,fixPar,betaCons,cons$deltaCons,hierStates,formula,formulaDelta,data,mixtures,field="fixPar")
    par  <- mapPar(hierBeta,hierDelta,fixPar,betaCons,cons$deltaCons,hierStates,formula,formulaDelta,data,mixtures,field="beta")
    wb   <- mapBounds(hierBeta,hierDelta,fixPar,betaCons,cons$deltaCons,hierStates,formula,formulaDelta,data,mixtures)
    
    betaCons <- cons$betaCons
    deltaCons <- cons$deltaCons
    
    fixPar$beta <- fix$beta
    fixPar$delta <- fix$delta
    
    beta0 <- par$beta
    delta0 <- par$delta
    
    workBounds <- list(beta=wb$beta, delta=wb$delta)
    if(is.null(workBounds$beta)){
      if(any(is.finite(betaLower))){
        workBounds$beta <- cbind(c(betaLower),c(betaUpper))
      }
    }
    
    if(!is.null(beta0)) {
      if(is.null(hierDelta)) beta0[which(is.na(beta0))] <- 0
    }
    
    if(mixtures>1) beta0 <- list(beta = beta0, pi = hierBeta$pi)
    else {
      if(!is.null(fixPar$beta)) rownames(fixPar$beta) <- colnames(covs)
      if(!is.null(fixPar$delta)) rownames(fixPar$delta) <- colnames(covsDelta)
      if(!is.null(betaCons)) rownames(betaCons) <- rownames(fixPar$beta)
      if(!is.null(deltaCons)) rownames(deltaCons) <- colnames(covsDelta)
      if(!is.null(beta0)) rownames(beta0) <- rownames(fixPar$beta)
      if(!is.null(delta0)) rownames(delta0) <- rownames(fixPar$delta)
    }
    
    # populate hierBeta and hierDelta if not provided
    hier <- mapHier(beta0,Pi,delta0,hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,g0,theta)
    hierBeta <- hier$hierBeta
    hierDelta <- hier$hierDelta
  } else betaCons <- deltaCons <- fixPar <- workBounds <- NULL
  
  return(list(nbStates=nbStates,dist=dist,formula=formula,formulaDelta=formulaDelta,beta=beta0,delta=delta0,hierBeta=hierBeta,hierDelta=hierDelta,betaRef=betaRef,betaCons=betaCons,deltaCons=deltaCons,fixPar=fixPar,workBounds=workBounds,stateNames=stateNames))
}

checkHierFormula <- function(data,hierFormula,hierStates,hierDist,checkData,what="formula"){
  if(is.null(hierFormula)){
    whierFormula <- data.tree::Node$new(hierStates$Get("name",filterFun=isRoot))
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
    if(!(what %in% hierFormula$fieldsAll)) stop(ifelse(what=="formula","'hierFormula'","'hierFormulaDelta'")," must include a ",ifelse(what=="formula","'formula'","'formulaDelta'")," field")
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

checkField <- function(what,field,level,hierBeta,hierStates,betaRef,nbCovs,mixtures,initial=FALSE,bounds=FALSE){
  if(level>1){
    #if(hierBeta[[paste0("level",level)]]$count!=hierStates$Get("count",filterFun = function(x) x$level==(level-1))) stop(what,"$",field," for level",level," must consist of ",hierStates$Get("count",filterFun = function(x) x$level==(level-1))," children: ",paste0(hierStates$Get("name",filterFun=function(x) x$level==level),collapse=", "))
    for(jj in hierStates$Get("name",filterFun=function(x) x$level==level)){
      inits <- hierBeta[[paste0("level",level)]][[jj]][[field]]
      nlStates <- hierStates$Get("count",filterFun=function(x) x$level==level)[jj]
      jRef <- sum(!(rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==level)[jj],times=hierStates$Get("leafCount",filterFun=function(x) x$level==level)[jj]) %in% betaRef))
      if(is.null(inits) & nlStates) stop(what,"$level",level,"$",jj,"$",field," are missing")
      pCount <- ifelse(initial,(nlStates-1),nlStates*(nlStates-1))
      if(!bounds) {
        if(nlStates && (!is.matrix(inits) || (ncol(inits)!=pCount*ifelse(jRef,jRef,1) | nrow(inits)!=nbCovs*mixtures))) stop(what,"$level",level,"$",jj,"$",field," should consist of ",nbCovs*mixtures," rows and ",pCount*ifelse(jRef,jRef,1), " columns")
      } else {
        if(nlStates && (!is.matrix(inits) || (nrow(inits)!=pCount*ifelse(jRef,jRef,1)*nbCovs*mixtures | ncol(inits)!=2))) stop(what,"$level",level,"$",jj,"$",field," should consist of ",pCount*ifelse(jRef,jRef,1)*nbCovs*mixtures," rows and 2 columns")
      }
    }
    inits <- hierBeta[[paste0("level",level)]]$Get(field,filterFun=isLeaf,simplify=FALSE)[hierStates$Get("name",filterFun=function(x) x$level==level)]
    if(!bounds) inits <- unlist(inits)
  } else {
    inits <- hierBeta[[paste0("level",level)]][[field]]
    if(is.null(inits)) stop(what,"$",field," are missing for level",level)
    pCount <- ifelse(initial,(hierStates$count-1),hierStates$count*(hierStates$count-1))
    if(!bounds) {
      if(!is.matrix(inits) || (ncol(inits)!=pCount | nrow(inits)!=nbCovs*mixtures)) stop(what,"$level",level,"$",field," should consist of ",nbCovs*mixtures," rows and ",pCount, " columns")
    } else {
      if(!is.matrix(inits) || (nrow(inits)!=pCount*nbCovs*mixtures | ncol(inits)!=2)) stop(what,"$level",level,"$",field," should consist of ",pCount*nbCovs*mixtures," rows and 2 columns")
    }
    inits <- hierBeta[[paste0("level",level)]][[field]]  
  }
  inits
}

mapCons <- function(hierBeta,hierDelta,fixPar,betaCons,hierStates,formula,formulaDelta,data,mixtures){
  
  bc <- betaCons
  deltaCons <- matrix(1:length(fixPar$delta),nrow(fixPar$delta),ncol(fixPar$delta),dimnames=dimnames(fixPar$delta))
  fixInd <- which(fixPar$delta==-1.e+10)
  if(length(fixInd)) deltaCons[fixInd] <- fixInd[1]
  
  betaRef <- rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  if("betaCons" %in% hierBeta$fieldsAll){
    what <- "hierBeta"
    field <- "betaCons"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)

      inits <- checkField(what,field,j,hierBeta,hierStates,betaRef,nbCovs,mixtures)
      betaInd <- betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]
      naInd <- is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),])
      bc[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][naInd] <- unique(betaInd[naInd])[inits][match(betaInd[naInd],unique(betaInd[naInd]))]

    }
  }
  if("deltaCons" %in% hierDelta$fieldsAll){
    what <- "hierDelta"
    field <- "deltaCons"
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        covNames <- colnames(model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
        covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE)
        betaInd <- betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]
        naInd <- is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),])
        bc[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][naInd] <- unique(betaInd[naInd])[inits][match(betaInd[naInd],unique(betaInd[naInd]))]
      } else if(j==1){
        covNames <- colnames(model.matrix(formulaDelta,data))
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE)
        initsInd <- which(is.na(fixPar$delta))#[(mix-1)*nbCovs+1:nbCovs,,drop=FALSE]))
        deltaCons[initsInd] <- initsInd[inits]
      }
    }
  }
  return(list(betaCons=bc, deltaCons = deltaCons))
}

mapPar <- function(hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,field="beta"){
  
  match.arg(field,c("beta","fixPar"))
  
  betaRef <- rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  if(field=="beta") beta <- delta <- NULL
  else {
    beta <- fixPar$beta
    delta <- fixPar$delta
  }
  if(field %in% hierBeta$fieldsAll){
    beta <- fixPar$beta
    what <- "hierBeta"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)
      
      inits <- checkField(what,field,j,hierBeta,hierStates,betaRef,nbCovs,mixtures)
      initsInd <- unique(betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][which(is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]))])
      beta[initsInd] <- inits
    }
  }
  field <- ifelse(field=="beta","delta",field)
  if(field %in% hierDelta$fieldsAll){
    if(is.null(beta)) beta <- fixPar$beta
    delta <- fixPar$delta
    what <- "hierDelta"
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        covNames <- colnames(model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
        covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE)
        initsInd <- unique(betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][which(is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]))])
        beta[initsInd] <- inits
      } else if(j==1){
        covNames <- colnames(model.matrix(formulaDelta,data))
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE)
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

mapBounds <- function(hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,field="workBounds"){
  
  match.arg(field,"workBounds")
  
  betaRef <- rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  beta <- delta <- NULL
  
  if(field %in% hierBeta$fieldsAll){
    betaLower <- matrix(-Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
    betaUpper <- matrix(Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
    what <- "hierBeta"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)
      
      inits <- checkField(what,field,j,hierBeta,hierStates,betaRef,nbCovs,mixtures,bounds=TRUE)
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

  if(field %in% hierDelta$fieldsAll){
    if(is.null(beta)) {
      betaLower <- matrix(-Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
      betaUpper <- matrix(Inf,nrow(fixPar$beta),ncol(fixPar$beta),dimnames = dimnames(fixPar$beta))
    }
    deltaLower <- matrix(-Inf,nrow(fixPar$delta),ncol(fixPar$delta),dimnames = dimnames(fixPar$delta))
    deltaUpper <- matrix(Inf,nrow(fixPar$delta),ncol(fixPar$delta),dimnames = dimnames(fixPar$delta))
    
    what <- "hierDelta"
    for(j in 1:(hierStates$height-1)){
      if(j>1){
        covNames <- colnames(model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
        covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
        nbCovs <- length(covNames)

        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,bounds=TRUE)
        initsInd <- unique(betaCons[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),][which(is.na(fixPar$beta[paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs)),]))])
        betaLower[initsInd] <- unlist(lapply(inits,function(x) x[,1]))
        betaUpper[initsInd] <- unlist(lapply(inits,function(x) x[,2]))
      } else if(j==1){
        covNames <- colnames(model.matrix(formulaDelta,data))
        nbCovs <- length(covNames)
        
        inits <- checkField(what,field,j,hierDelta,hierStates,betaRef,nbCovs,mixtures,initial=TRUE,bounds=TRUE)
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

mapHier <- function(beta,pi,delta,hierBeta,hierDelta,fixPar,betaCons,deltaCons,hierStates,formula,formulaDelta,data,mixtures,g0=NULL,theta=NULL,fill=FALSE){
  
  if(is.null(hierBeta)){
    hierBeta <- Node$new("hierBeta")
    hierBeta$AddChild(paste0("level1"))
    for(j in 2:(hierStates$height-1)){
      hierBeta$AddChild(paste0("level",j))
      for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
        hierBeta[[paste0("level",j)]]$AddChild(jj)
      }
    }
  }
  if(is.null(hierDelta)){
    hierDelta <- Node$new("hierDelta")
    hierDelta$AddChild(paste0("level1"))
    for(j in 2:(hierStates$height-1)){
      hierDelta$AddChild(paste0("level",j))
      for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
        hierDelta[[paste0("level",j)]]$AddChild(jj)
      }
    }        
  }
  
  if(is.null(beta)) beta <- matrix(0,nrow(fixPar$beta),ncol(fixPar$beta),dimnames=dimnames(fixPar$beta))
  if(is.null(delta)) delta <- matrix(0,nrow(fixPar$delta),ncol(fixPar$delta),dimnames=dimnames(fixPar$delta))
  
  betaRef <- rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  
  whierBeta <- NULL
  if(is.list(hierBeta)){
    if(inherits(hierBeta$beta,"Node")) whierBeta <- Clone(hierBeta$beta)
  } else if(inherits(hierBeta,"Node")) whierBeta <- data.tree::Clone(hierBeta)

  if(inherits(whierBeta,"Node")){
    what <- "hierBeta"
    for(j in 1:(hierStates$height-1)){
      covNames <- colnames(model.matrix(formula,data[which(data$level==j),]))
      covNames <- covNames[grepl(paste0("level",j,"$"),covNames) | grepl(paste0("I((level == \"",j,"\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)
      if(mixtures>1) covNames <- paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs))
      if(j>1){
        initsInd <- unique(betaCons[covNames,][which(is.na(fixPar$beta[covNames,]))])
        dimNames <- mapply(function(x) mapply(`[[`, dimnames(fixPar$beta), arrayInd(x, dim(fixPar$beta))),initsInd)
        inits <- beta[initsInd]
        count <- 0
        t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
        names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
        for(jj in names(t)){
          nStates <- t[[jj]]$count#length(t[[jj]]$Get("state",filterFun = data.tree::isLeaf))
          jRef <- sum(!(rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j)[jj],times=hierStates$Get("leafCount",filterFun=function(x) x$level==j)[jj]) %in% betaRef))
          pCount <- nStates*(nStates-1)*ifelse(jRef,jRef,1)
          if(nStates){
            if(is.null(whierBeta[[paste0("level",j)]][[jj]]$betaCons)) {
              iRef <- count + matrix(1:(pCount*nbCovs*mixtures),nbCovs*mixtures,pCount)
            } else iRef <- whierBeta[[paste0("level",j)]][[jj]]$betaCons
            idimNames <- list(covNames,unique(dimNames[2,iRef]))
            if(fill & is.null(whierBeta[[paste0("level",j)]][[jj]]$betaCons)) whierBeta[[paste0("level",j)]][[jj]]$betaCons <- matrix(iRef,nbCovs*mixtures,pCount,dimnames=idimNames) 
            if(fill & is.null(whierBeta[[paste0("level",j)]][[jj]]$fixPar)) whierBeta[[paste0("level",j)]][[jj]]$fixPar <- matrix(NA,nbCovs*mixtures,pCount,dimnames=idimNames) 
            if(fill & is.null(whierBeta[[paste0("level",j)]][[jj]]$workBounds)) whierBeta[[paste0("level",j)]][[jj]]$workBounds <- matrix(c(-Inf,Inf),nbCovs*mixtures*pCount,2,byrow=TRUE, dimnames = list(paste0(idimNames[[1]],":",rep(idimNames[[2]],each=length(idimNames[[1]]))),c("lower","upper")))
            whierBeta[[paste0("level",j)]][[jj]]$beta <- matrix(inits[iRef],nbCovs*mixtures,pCount,dimnames=idimNames)
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
        if(fill & is.null(whierBeta[[paste0("level",j)]]$betaCons)) whierBeta[[paste0("level",j)]]$betaCons <- matrix(1:length(beta[initsInd]),nbCovs*mixtures,dimnames = dimNames)
        if(fill & is.null(whierBeta[[paste0("level",j)]]$fixPar)) whierBeta[[paste0("level",j)]]$fixPar <- matrix(rep(NA,length(beta[initsInd])),nbCovs*mixtures,dimnames = dimNames)
        if(fill & is.null(whierBeta[[paste0("level",j)]]$workBounds)) whierBeta[[paste0("level",j)]]$workBounds <- matrix(c(-Inf,Inf),length(beta[initsInd]),2,byrow=TRUE, dimnames = list(paste0(dimNames[[1]],":",rep(dimNames[[2]],each=length(dimNames[[1]]))),c("lower","upper")))
        whierBeta[[paste0("level",j)]]$beta <- matrix(beta[initsInd],nbCovs*mixtures,dimnames = dimNames)
      }
    }
  }
  
  whierDelta <- data.tree::Clone(hierDelta)
  what <- "hierDelta"
  for(j in 1:(hierStates$height-1)){
    if(j>1){
      covNames <- colnames(model.matrix(formula,data[which(data$level==paste0(j,"i")),]))
      covNames <- covNames[grepl(paste0("level",j,"i$"),covNames) | grepl(paste0("I((level == \"",j,"i\")"),covNames,fixed=TRUE)]
      nbCovs <- length(covNames)
      if(mixtures>1) covNames <- paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs))
      initsInd <- unique(betaCons[covNames,][which(is.na(fixPar$beta[covNames,]))])
      dimNames <- mapply(function(x) mapply(`[[`, dimnames(fixPar$beta), arrayInd(x, dim(fixPar$beta))),initsInd)
      inits <- beta[initsInd]
      count <- 0
      t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
      names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
      for(jj in names(t)){
        nStates <- nStates <- t[[jj]]$count#length(t[[jj]]$Get("state",filterFun = data.tree::isLeaf))
        jRef <- sum(!(rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j)[jj],times=hierStates$Get("leafCount",filterFun=function(x) x$level==j)[jj]) %in% betaRef))
        pCount <- (nStates-1)*ifelse(jRef,jRef,1)
        if(nStates){
          if(is.null(whierDelta[[paste0("level",j)]][[jj]]$deltaCons)) {
            iRef <- count + matrix(1:(pCount*nbCovs*mixtures),nbCovs*mixtures,pCount)
          } else iRef <- whierDelta[[paste0("level",j)]][[jj]]$deltaCons
          idimNames <- list(covNames,gsub('.*->',"state",unique(dimNames[2,iRef])))
          if(fill & is.null(whierDelta[[paste0("level",j)]][[jj]]$deltaCons)) whierDelta[[paste0("level",j)]][[jj]]$deltaCons <- matrix(iRef,nbCovs*mixtures,pCount,dimnames=idimNames)
          if(fill & is.null(whierDelta[[paste0("level",j)]][[jj]]$fixPar)) whierDelta[[paste0("level",j)]][[jj]]$fixPar <- matrix(NA,nbCovs*mixtures,pCount,dimnames=idimNames)
          if(fill & is.null(whierDelta[[paste0("level",j)]][[jj]]$workBounds)) whierDelta[[paste0("level",j)]][[jj]]$workBounds <- matrix(c(-Inf,Inf),nbCovs*mixtures*pCount,2,byrow=TRUE, dimnames = list(paste0(idimNames[[1]],":",rep(idimNames[[2]],each=length(idimNames[[1]]))),c("lower","upper")))
          whierDelta[[paste0("level",j)]][[jj]]$delta <- matrix(inits[iRef],nbCovs*mixtures,pCount,dimnames=idimNames)
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
      covNames <- colnames(model.matrix(formulaDelta,data))
      nbCovs <- length(covNames)
      if(mixtures>1) covNames <- paste0(covNames,"_mix",rep(1:mixtures,each=nbCovs))
      initsInd <- unique(deltaCons[covNames,][which(is.na(fixPar$delta[covNames,]))])
      dimNames <- mapply(function(x) mapply(`[[`, dimnames(fixPar$delta), arrayInd(x, dim(fixPar$delta))),initsInd)
      dimNames <- list(covNames,unique(dimNames[2,]))
      if(fill & is.null(whierDelta[[paste0("level",j)]]$deltaCons)) whierDelta[[paste0("level",j)]]$deltaCons <- matrix(1:length(delta[initsInd]),nbCovs*mixtures,dimnames=dimNames)
      if(fill & is.null(whierDelta[[paste0("level",j)]]$fixPar)) whierDelta[[paste0("level",j)]]$fixPar <- matrix(rep(NA,length(delta[initsInd])),nbCovs*mixtures,dimnames=dimNames)
      if(fill & is.null(whierDelta[[paste0("level",j)]]$workBounds)) whierDelta[[paste0("level",j)]]$workBounds <- matrix(c(-Inf,Inf),length(delta[initsInd]),2,byrow=TRUE, dimnames = list(paste0(dimNames[[1]],":",rep(dimNames[[2]],each=length(dimNames[[1]]))),c("lower","upper")))
      whierDelta[[paste0("level",j)]]$delta <- matrix(delta[initsInd],nbCovs*mixtures,dimnames=dimNames)
    }
  }
  
  recharge <- newFormulas(formula,length(hierStates$Get("state",filterFun=data.tree::isLeaf)))$recharge
  if(mixtures>1 | !is.null(recharge)){
    whierBeta <- list(beta=whierBeta)
    if(mixtures>1) whierBeta$pi <- pi
    if(!is.null(recharge)){
      whierBeta$g0 <- g0
      whierBeta$theta <- theta
    }
  }
  return(list(hierBeta=whierBeta,hierDelta=whierDelta))
}
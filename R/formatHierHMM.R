#' @importFrom data.tree Node Get Do ToDataFrameTypeCol Traverse Aggregate AreNamesUnique isRoot isLeaf
#' @importFrom stats terms
formatHierHMM <- function(data=NULL,hierStates,hierDist,
                          hierFormula,formulaDelta,mixtures,
                          workBounds,betaCons,
                          fixPar,checkData=TRUE){
  
  if(!inherits(hierStates,"Node")) stop("'hierStates' must be of class Node; see ?data.tree::Node")
  if(!("state" %in% hierStates$fieldsAll)) stop("'hierStates' must include a 'state' field")
  nbLevels <- 2*hierStates$height - 3 #hierStates$height #ncol(data.tree::ToDataFrameTypeCol(hierStates))
  
  hdf <- data.tree::ToDataFrameTypeCol(hierStates, "state")
  if(any(is.na(hdf$state))) stop("'state' field in 'hierStates' cannot contain NAs")
  nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))
  
  if(any(duplicated(hdf$state))) stop("'state' field in 'hierStates' cannot contain duplicates")
  if(any(sort(hdf$state)!=1:nbStates)) stop("'state' field in 'hierStates' must include all integers between 1 and ",nbStates)
  if(!data.tree::AreNamesUnique(hierStates)) stop("node names in 'hierStates' must be unique")
  #if(any(is.na(hdf))) stop("missing levels are not permitted in 'hierStates'")
  
  stateNames <- unname(hierStates$Get("name",filterFun=data.tree::isLeaf)[hdf$state])
  
  if(hierStates$height<=2) stop("'hierStates' must contain >2 levels")
  
  if(!inherits(hierDist,"Node")) stop("'hierDist' must be of class Node; see ?data.tree::Node")
  if(!("dist" %in% hierDist$fieldsAll)) stop("'hierDist' must include a 'dist' field")
  if(!data.tree::AreNamesUnique(hierDist)) stop("node names in 'hierDist' must be unique")
  if(hierDist$height!=3) stop("'hierDist' hierarchy must contain 2 levels")
  
  if(!is.null(data) && checkData){
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
  
  if(is.null(hierFormula)){
    hierFormula <- data.tree::Node$new(hierStates$Get("name",filterFun=isRoot))
    hierFormula$AddChild(hierDist$Get("name",filterFun=function(x) x$level==2)[1],formula=~1)
    for(j in hierDist$Get("name",filterFun=function(x) x$level==2)[-1]){
      hierFormula$AddChild(paste0(j,"i"),formula=~1)
      hierFormula$AddChild(j,formula=~1)
    }
  } else {
    if(!inherits(hierFormula,"Node")) stop("'hierFormula' must be of class Node; see ?data.tree::Node")
    if(!("formula" %in% hierFormula$fieldsAll)) stop("'hierFormula' must include a 'formula' field")
    if(!data.tree::AreNamesUnique(hierFormula)) stop("node names in 'hierFormula' must be unique")
    if(hierFormula$height!=2) stop("'hierFormula' hierarchy must contain 1 level")
    
    if(!is.null(data) && checkData){
      if(!all(hierFormula$Get("name",filterFun=function(x) x$level==2) %in% paste0("level",levels(data$level)))) 
        stop("'hierFormula' level types can only include ",paste(paste0("level",levels(data$level)),collapse=", "))
    }
    
    #if(!all(hierDist$Get("name",filterFun=function(x) x$level==2)==hierFormula$Get("name",filterFun=function(x) x$level==2)[seq(1,hierFormula$count,2)])) stop("'hierDist' and 'hierFormula' are not consistent; check number of nodes, node names, and node order")
    for(j in hierDist$Get("name",filterFun=function(x) x$level==2)){
      if(is.null(hierFormula[[j]])) {
        hierFormula$AddChild(j,formula=~1)
      } else if(is.null(hierFormula[[j]]$formula)){
        hierFormula[[j]]$formula <- ~1
      } else if(!is.formula(hierFormula[[j]]$formula)) stop("'hierFormula$",j,"$formula' must be a formula")
      
      if(j!=hierDist$Get("name",filterFun=function(x) x$level==2)[1]){
        if(is.null(hierFormula[[paste0(j,"i")]])) {
          hierFormula$AddChild(paste0(j,"i"),formula=~1)
        } else if(is.null(hierFormula[[paste0(j,"i")]]$formula)){
          hierFormula[[paste0(j,"i")]]$formula <- ~1
        } else if(!is.formula(hierFormula[[paste0(j,"i")]]$formula)) stop("'hierFormula$",paste0(j,"i"),"$formula' must be a formula")
      }
    }
  }
  if(is.null(formulaDelta)) stop("formulaDelta cannot be NULL")
  else {
    deltaTerms <- attr(stats::terms(formulaDelta),"term.labels")
    if(length(deltaTerms)){
      if(grepl("level",deltaTerms)) stop("'level' cannot be included in formulaDelta")
    }
    if(!attr(stats::terms(formulaDelta),"intercept")) stop("formulaDelta must include an intercept term")
  }
  #if(attr(stats::terms(formula),"intercept")) stop("formula can not include an intercept term")

  formula <- formatHierFormula(hierFormula)
  
  # set t.p.m. reference states based on top level
  #if(is.null(betaRef)){
    betaRef <- rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  #}
  
  if(!is.null(data) && (is.null(fixPar$beta) | is.null(betaCons) | is.null(fixPar$delta))){
    
    covs <- model.matrix(formula,data)
    nbCovs <- ncol(covs)
    
    if(is.null(fixPar$beta)){
      betaInd <- 1
      if(is.null(fixPar)) fixPar <- list()
      fixPar$beta <- matrix(NA,nbCovs*mixtures,nbStates*(nbStates-1))
    } else betaInd <- 0
    if(!all(dim(fixPar$beta)==c(nbCovs*mixtures,nbStates*(nbStates-1)))) stop("fixPar$beta must be a matrix of dimension ",nbCovs*mixtures,"x",nbStates*(nbStates-1))
    rownames(fixPar$beta) <- paste0(rep(colnames(covs),mixtures),"_mix",rep(1:mixtures,each=length(colnames(covs))))
    colnames(fixPar$beta) <- c(sapply(1:nbStates,function(x) paste(rep(x,each=nbStates-1),"->",1:nbStates)[-betaRef[x]]))
    
    if(is.null(betaCons)){
      conInd <- 1
      betaCons <- matrix(1:(nbCovs*mixtures*nbStates*(nbStates-1)),nbCovs*mixtures,nbStates*(nbStates-1))
    } else conInd <- 0
    if(!all(dim(betaCons)==c(nbCovs*mixtures,nbStates*(nbStates-1)))) stop("betaCons must be a matrix of dimension ",nbCovs*mixtures,"x",nbStates*(nbStates-1))
    dimnames(betaCons)<-dimnames(fixPar$beta)
    
    if(is.null(fixPar$delta)){
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
    }
    
    hierStates$Do(function(node) node$levelStates <- node$Get("state",traversal="level",filterFun=data.tree::isLeaf))
    
    #top level t.p.m. (must map all states at bottom level back to top level)
    j <- 1
    t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
    names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
    levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-1],"\")"),colnames(covs),fixed=TRUE)]
    otherCovs <- colnames(covs)[-which(colnames(covs) %in% c(paste0("level",levels(data$level)),levelCovs))]
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
              if(conInd){
                betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),withinConstr[[h]]] <- min(betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),withinConstr[[h]]])
                for(lcovs in levelCovs){
                  betaCons[paste0(lcovs,"_mix",mix),withinConstr[[h]]] <- min(betaCons[paste0(lcovs,"_mix",mix),withinConstr[[h]]])
                }
              }
              if(betaInd){
                fixPar$beta[paste0("level",levels(data$level)[2*j-1],"_mix",mix),withinConstr[[h]]] <- -1.e+10
                for(lcovs in levelCovs){
                  fixPar$beta[paste0(lcovs,"_mix",mix),withinConstr[[h]]] <- -1.e+10
                }
              }
            }
          }
          # constrain across transitions to reference states
          levelStates <- unlist(lapply(tt[(1:length(tt))[-h]],function(x) x$Get("state",filterFun = data.tree::isLeaf)))
          acrossConstr[[h]] <- match(paste0(rep(levelStates,each=length(stateInd[-1]))," -> ",stateInd[-1]),colnames(fixPar$beta),nomatch=0)
          acrossRef <- match(paste0(levelStates," -> ",stateInd[1]),colnames(fixPar$beta),nomatch=0)
          for(mix in 1:mixtures){
            if(conInd){
              if(any(acrossConstr[[h]])){
                betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),acrossConstr[[h]]] <- min(betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),acrossConstr[[h]]])
                for(lcovs in levelCovs){
                  betaCons[paste0(lcovs,"_mix",mix),acrossConstr[[h]]] <- min(betaCons[paste0(lcovs,"_mix",mix),acrossConstr[[h]]])
                }
                betaCons[paste0("level",levels(data$level)[(2*j):nbLevels],"_mix",mix),acrossConstr[[h]]] <- min(betaCons[paste0("level",levels(data$level)[(2*j):nbLevels],"_mix",mix),acrossConstr[[h]]])
                for(ocovs in otherCovs){
                  betaCons[ocovs,acrossConstr[[h]]] <- min(betaCons[ocovs,acrossConstr[[h]]])
                }
              }
              if(any(acrossRef)){
                betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),acrossRef] <- min(betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),acrossRef])
                for(lcovs in levelCovs){
                  betaCons[paste0(lcovs,"_mix",mix),acrossRef] <- min(betaCons[paste0(lcovs,"_mix",mix),acrossRef])
                }
                betaCons[paste0("level",levels(data$level)[(2*j):nbLevels],"_mix",mix),acrossRef] <- min(betaCons[paste0("level",levels(data$level)[(2*j):nbLevels],"_mix",mix),acrossRef])
                for(ocovs in otherCovs){
                  betaCons[ocovs,acrossRef] <- min(betaCons[ocovs,acrossRef])
                }
              }
            }
            if(betaInd){
              if(any(acrossConstr[[h]])){
                fixPar$beta[paste0("level",levels(data$level)[2*j-1],"_mix",mix),acrossConstr[[h]]] <- -1.e+10
                for(lcovs in levelCovs){
                  fixPar$beta[paste0(lcovs,"_mix",mix),acrossConstr[[h]]] <- -1.e+10
                }
                fixPar$beta[paste0("level",levels(data$level)[(2*j):nbLevels],"_mix",mix),acrossConstr[[h]]] <- -1.e+10
                for(ocovs in otherCovs){
                  fixPar$beta[ocovs,acrossConstr[[h]]] <- -1.e+10
                }
              }
              if(any(acrossRef)){
                fixPar$beta[paste0("level",levels(data$level)[(2*j):nbLevels],"_mix",mix),acrossRef] <- -1.e+10
                for(ocovs in otherCovs){
                  fixPar$beta[ocovs,acrossRef] <- -1.e+10
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
      
      levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-2],"\")"),colnames(covs),fixed=TRUE)]
      otherCovs <- colnames(covs)[-which(colnames(covs) %in% c(paste0("level",levels(data$level)),levelCovs))]
      
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
            if(conInd){
              betaCons[paste0("level",levels(data$level)[2*j-2],"_mix",mix),initConstr] <- min(betaCons[paste0("level",levels(data$level)[2*j-2],"_mix",mix),initConstr])
              for(lcovs in levelCovs){
                betaCons[paste0(lcovs,"_mix",mix),initConstr] <- min(betaCons[paste0(lcovs,"_mix",mix),initConstr])
              }
            }
            if(betaInd){
              fixPar$beta[paste0("level",levels(data$level)[2*j-2],"_mix",mix),initConstr] <- -1.e+10
              for(lcovs in levelCovs){
                fixPar$beta[paste0(lcovs,"_mix",mix),initConstr] <- -1.e+10
              }
            }
          }
          if(!any(betaRef[fromState] %in% toStates)){
            #warning("to state = ",toStates,"; betaRef[fromState] = ",betaRef[fromState])
            if(any(refConstr)){
              if(sum(refConstr>0)>1){
                betaLower[paste0("level",levels(data$level)[2*j-2],"_mix",mix),refConstr] <- 500
                for(lcovs in levelCovs){
                  betaLower[paste0(lcovs,"_mix",mix),refConstr] <- -(500-100)/length(levelCovs) # don't let XB be less than 100
                }
              } else if(betaInd) {
                if(!length(levelCovs)) fixPar$beta[paste0("level",levels(data$level)[2*j-2],"_mix",mix),refConstr] <- 500
              }
            }
          }
        }
      }
      
      levelCovs <- colnames(covs)[grepl(paste0("I((level == \"",levels(data$level)[2*j-1],"\")"),colnames(covs),fixed=TRUE)]
      
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
            if(conInd){
              betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),initConstr] <- min(betaCons[paste0("level",levels(data$level)[2*j-1],"_mix",mix),initConstr])
              for(lcovs in levelCovs){
                betaCons[paste0(lcovs,"_mix",mix),initConstr] <- min(betaCons[paste0(lcovs,"_mix",mix),initConstr])
              }
            }
            if(betaInd){
              fixPar$beta[paste0("level",levels(data$level)[2*j-1],"_mix",mix),initConstr] <- -1.e+10
              for(lcovs in levelCovs){
                fixPar$beta[paste0(lcovs,"_mix",mix),initConstr] <- -1.e+10
              }
            }
          }
          if(!any(betaRef[fromStates] %in% toStates)){
            #warning("to state = ",toStates,"; betaRef[fromStates] = ",betaRef[fromStates])
            if(any(refConstr)) {
              if(sum(refConstr>0)>1){
                betaLower[paste0("level",levels(data$level)[2*j-1],"_mix",mix),refConstr] <- 500
                for(lcovs in levelCovs){
                  betaLower[paste0(lcovs,"_mix",mix),refConstr] <- -(500-100)/length(levelCovs) # don't let XB be less than 100
                }
              } else if(betaInd) {
                if(!length(levelCovs)) fixPar$beta[paste0("level",levels(data$level)[2*j-1],"_mix",mix),refConstr] <- 500
              }
            }
          }
        }
      }
    }
    if(conInd){
      fixInd <- which(fixPar$beta==-1.e+10)
      if(length(fixInd)) betaCons[fixInd] <- fixInd[1]
    }
    if(is.null(workBounds$beta) & conInd){
      if(any(is.finite(betaLower))){
        workBoundsBeta <- cbind(c(betaLower),c(betaUpper))
        if(is.null(workBounds)) workBounds <- list()
        workBounds$beta <- workBoundsBeta
      }
    }
    if(mixtures==1){
      rownames(fixPar$delta) <- colnames(covsDelta)
      rownames(fixPar$beta) <- colnames(covs)
      rownames(betaCons) <- rownames(fixPar$beta)
    }
    return(list(nbStates=nbStates,dist=dist,formula=formula,betaRef=betaRef,betaCons=betaCons,fixPar=fixPar,workBounds=workBounds,stateNames=stateNames))
  } else {
    return(list(nbStates=nbStates,dist=dist,formula=formula,betaRef=betaRef,stateNames=stateNames))
  }
}
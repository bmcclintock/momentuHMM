#' @importFrom stats terms
#' @importFrom data.tree Get
formatHierFormula <- function(data,hierFormula,hierStates){
  lLevels <- sort(hierFormula$Get("name",filterFun=function(x) x$level==2))[c(1,1+rep(2:1,hierStates$height-2)+rep(seq(0,(hierStates$height-3)*2,2),each=2))]
  formulaTerms <- formTerms <- recTerms <- list()
  recharge <- NULL
  betaRef <- rep(hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2),times=hierStates$Get("leafCount",filterFun=function(x) x$level==2))
  if(!is.null(data)) factorTerms <- names(data)[which(unlist(lapply(data,function(x) inherits(x,"factor"))))]
  for(j in lLevels){
    newForm <- newFormulas(hierFormula[[j]]$formula, length(hierStates$Get("state",filterFun=data.tree::isLeaf)), betaRef, hierarchical=TRUE)
    formulaTerms[[j]] <- newForm$formterms
    if(any(grepl("level",formulaTerms[[j]]))) stop("hierFormula$",j,"$formula cannot include 'level'")
    if(is.null(data)) {
      formTerms[[j]] <- paste0("I((level=='",gsub("level","",j),"')*1)")
      if(length(formulaTerms[[j]])) formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*",formulaTerms[[j]],")"))
    } else {
      formTerms[[j]] <- getFactorTerms(formulaTerms[[j]],factorTerms,newForm$newformula,data,j)
    }
    if(!is.null(newForm$recharge)){
      recharge[[j]] <- newForm$recharge
      for(parm in c("g0","theta")){
        rechargeTerms <- attr(stats::terms(recharge[[j]][[parm]]),"term.labels")
        if(length(rechargeTerms)){
          if(any(grepl("level",rechargeTerms))) stop("hierFormula$",j," recharge formula cannot include 'level'")
          if(is.null(data)) recTerms <- paste0("I((level=='",gsub("level","",j),"')*1):",as.character(recharge[[j]][[parm]])[-1])
          else recTerms <- getFactorTerms(rechargeTerms,factorTerms,recharge[[j]][[parm]],data,j)
        } else {
          if(parm=="theta" & !attributes(stats::terms(recharge[[j]]$theta))$intercept) stop("invalid recharge model for ",j," -- theta must include an intercept and at least 1 covariate")
          recTerms <- paste0("I((level=='",gsub("level","",j),"')*1)")
        }
        form <- paste0("~0+",paste0(unlist(recTerms),collapse = " + "))
        recharge[[j]][[parm]] <- stats::as.formula(form)
      }
      formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*1):recharge(g0=~",as.character(recharge[[j]]$g0)[-1],", theta=~",as.character(recharge[[j]]$theta)[-1],")"))
      if(!is.null(data)) data[[paste0("recharge",gsub("level","",j))]] <- rep(0,nrow(data))
    }
  }
  form <- paste0("~ 0 + ",paste0(unlist(formTerms),collapse = " + "))
  form <- stats::as.formula(form)
  return(list(formula=form,data=data,recharge=recharge))
}

getFactorTerms <- function(formulaTerms,factorTerms,formula,data,level){
  mm<-stats::model.matrix(formula,data)
  as <- attr(mm,"assign")
  colmm <- colnames(mm)
  for(jj in 1:length(colmm)){
    cterms <- unlist(strsplit(formulaTerms[as[jj]],":"))
    cInd <- which(cterms %in% factorTerms)
    lInd <- unlist(strsplit(colmm[jj],":"))[cInd]
    if(length(cInd)){
      for(kk in 1:length(cInd)){
        lInd[kk] <- gsub(cterms[cInd[kk]],"",lInd[kk])
      }
    }
    names(lInd) <- cterms[cInd]
    for(kk in cterms[cInd]){
      colmm[jj] <- gsub(paste0(kk,lInd[kk]),paste0("(",kk,"=='",lInd[kk],"')"),colmm[jj])
    }
    colmm[jj] <- paste0("I((level=='",gsub("level","",level),"')*",colmm[jj],")")
    colmm[jj] <- gsub(":","*",colmm[jj])
    colmm[jj] <- gsub("*(Intercept)","*1",colmm[jj],fixed=TRUE)
  }
  colmm
}

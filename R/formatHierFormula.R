#' @importFrom stats terms
#' @importFrom data.tree Get
formatHierFormula <- function(data,hierFormula,hierStates){
  if(is.null(data)){
    lLevels <- hierFormula$Get("name",filterFun=function(x) x$level==2)
    formulaTerms <- list()
    recharge <- NULL
    for(j in lLevels){
      newForm <- newFormulas(hierFormula[[j]]$formula, length(hierStates$Get("state",filterFun=data.tree::isLeaf)),hierarchical=TRUE)
      formulaTerms[[j]] <- newForm$formterms
      if(length(formulaTerms[[j]])){
        if(any(grepl("level",formulaTerms[[j]]))) stop("hierFormula$",j,"$formula cannot include 'level'")
        formulaTerms[[j]] <- paste0("I((level=='",gsub("level","",j),"')*",formulaTerms[[j]],")")
      }
      if(!is.null(newForm$recharge)){
        recharge[[j]] <- newForm$recharge
        for(parm in c("g0","theta")){
          rechargeTerms <- attr(stats::terms(recharge[[j]][[parm]]),"term.labels")
          if(length(rechargeTerms)){
            if(any(grepl("level",rechargeTerms))) stop("hierFormula$",j," recharge formula cannot include 'level'")
            recTerms <- paste0("I((level=='",gsub("level","",j),"')*1):",as.character(recharge[[j]][[parm]])[-1])
          } else {
            if(parm=="theta" & !attributes(terms(recharge[[j]]$theta))$intercept) stop("invalid recharge model for ",j," -- theta must include an intercept and at least 1 covariate")
            recTerms <- paste0("I((level=='",gsub("level","",j),"')*1)")
          }
          form <- paste0("~0+",paste0(unlist(recTerms),collapse = " + "))
          recharge[[j]][[parm]] <- as.formula(form)
        }
        formulaTerms[[j]] <- c(formulaTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*1):recharge(g0=~",as.character(recharge[[j]]$g0)[-1],", theta=~",as.character(recharge[[j]]$theta)[-1],")"))
      }
    }
    form <- "~ 0 + level"
    if(length(unlist(formulaTerms))) form <- paste0(form," + ",paste0(unlist(formulaTerms),collapse = " + "))
    form <- as.formula(form)
  } else {
    lLevels <- paste0("level",levels(data$level))#hierFormula$Get("name",filterFun=function(x) x$level==2)
    formulaTerms <- formTerms <- recTerms <- list()
    recharge <- NULL
    factorTerms <- names(data)[which(unlist(lapply(data,function(x) inherits(x,"factor"))))]
    for(j in lLevels){
      newForm <- newFormulas(hierFormula[[j]]$formula, length(hierStates$Get("state",filterFun=data.tree::isLeaf)),hierarchical=TRUE)
      #if(j==lLevels[1]) if(!attr(stats::terms(hierFormula[[j]]$formula),"intercept")) stop("hierFormula$",j,"$formula must include an intercept term")
      formulaTerms[[j]] <- newForm$formterms
      if(length(formulaTerms[[j]])){
        if(any(grepl("level",formulaTerms[[j]]))) stop("hierFormula$",j,"$formula cannot include 'level'")
        formTerms[[j]] <- getFactorTerms(formulaTerms[[j]],factorTerms,hierFormula[[j]]$formula,data,j)
      }
      if(!is.null(newForm$recharge)){
        recharge[[j]] <- newForm$recharge
        for(parm in c("g0","theta")){
          rechargeTerms <- attr(stats::terms(recharge[[j]][[parm]]),"term.labels")
          if(length(rechargeTerms)){
            if(any(grepl("level",rechargeTerms))) stop("hierFormula$",j," recharge formula cannot include 'level'")
            recTerms <- getFactorTerms(rechargeTerms,factorTerms,recharge[[j]][[parm]],data,j)
          } else {
            if(parm=="theta" & !attributes(terms(recharge[[j]]$theta))$intercept) stop("invalid recharge model for ",j," -- theta must include an intercept and at least 1 covariate")
            recTerms <- paste0("I((level=='",gsub("level","",j),"')*1)")
          }
          form <- paste0("~0+",paste0(unlist(recTerms),collapse = " + "))
          recharge[[j]][[parm]] <- as.formula(form)
        }
        formulaTerms[[j]] <- c(formulaTerms[[j]],"recharge")
        formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*1):recharge(g0=~",as.character(recharge[[j]]$g0)[-1],", theta=~",as.character(recharge[[j]]$theta)[-1],")"))
        data$recharge <- rep(0,nrow(data))
      }
    }
    factorLevels <- which(unlist(lapply(formulaTerms,function(x) any(x %in% factorTerms))))
    if(any(factorLevels)){
      if(length(factorLevels)==length(lLevels)) {
        form <- "~ 0"
      } else form <- paste0("~ 0 + ",paste0("I((level=='",gsub("level","",lLevels[-factorLevels]),"')*1)",collapse=" + "))
      intercepts <- paste0("I((level=='",gsub("level","",lLevels),"')*1)")
    } else {
      form <- paste0("~ 0 + ",paste0("I((level=='",gsub("level","",lLevels),"')*1)",collapse=" + "))
      intercepts <- paste0("I((level=='",gsub("level","",lLevels),"')*1)")
      #form <- "~ 0 + level"
      #intercepts <- paste0("level",levels(data$level))
    }
    if(length(unlist(formulaTerms))) form <- paste0(form," + ",paste0(unlist(formTerms),collapse = " + "))
    form <- as.formula(form)
  }
  return(list(formula=form,data=data,recharge=recharge))
}

getFactorTerms <- function(formulaTerms,factorTerms,formula,data,level){
  if(any(formulaTerms %in% factorTerms)){
    mm<-model.matrix(formula,data)
    as <- attr(mm,"assign")
    colmm <- colnames(mm)
    for(jj in 1:length(colmm)){
      cterms <- unlist(strsplit(formulaTerms[as[jj]],":"))
      lInd <- unlist(strsplit(colmm[jj],":"))[which(cterms %in% factorTerms)]
      for(kk in cterms[which(cterms %in% factorTerms)]){
        lInd <- gsub(kk,"",lInd)
      }
      names(lInd) <- cterms[which(cterms %in% factorTerms)]
      for(kk in cterms[which(cterms %in% factorTerms)]){
        colmm[jj] <- gsub(paste0(kk,lInd[kk]),paste0("(",kk,"=='",lInd[kk],"')"),colmm[jj])
      }
      colmm[jj] <- paste0("I((level=='",gsub("level","",level),"')*",colmm[jj],")")
      colmm[jj] <- gsub(":","*",colmm[jj])
      colmm[jj] <- gsub("*(Intercept)","*1",colmm[jj],fixed=TRUE)
    }
    formTerms <- colmm
  } else {
    for(k in 1:length(formulaTerms)){
      formTerms <- paste0("I((level=='",gsub("level","",level),"')*",formulaTerms,")")
    }
  }
  formTerms
}
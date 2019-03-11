#' @importFrom stats terms
#' @importFrom data.tree Get
formatHierFormula <- function(data,hierFormula,hierStates){
  if(is.null(data)){
    lLevels <- hierFormula$Get("name",filterFun=function(x) x$level==2)
    formulaTerms <- recharge <- list()
    for(j in lLevels){
      newForm <- newFormulas(hierFormula[[j]]$formula, length(hierStates$Get("state",filterFun=data.tree::isLeaf)))
      formulaTerms[[j]] <- newForm$formterms
      recharge[[j]] <- newForm$recharge
      if(length(formulaTerms[[j]])){
        if(any(grepl("level",formulaTerms[[j]]))) stop("hierFormula$",j,"$formula cannot include 'level'")
        formulaTerms[[j]] <- paste0("I((level=='",gsub("level","",j),"')*",formulaTerms[[j]],")")
      }
      if(!is.null(recharge[[j]])){
        #formulaTerms[[j]] <- c(formulaTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*","recharge",gsub("level","",j),")"))
        #data[[paste0("recharge",gsub("level","",j))]] <- rep(0,nrow(data))
        formulaTerms[[j]] <- c(formulaTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*","recharge)"))
        data$recharge <- rep(0,nrow(data))
      }
    }
    form <- "~ 0 + level"
    if(length(unlist(formulaTerms))) form <- paste0(form," + ",paste0(unlist(formulaTerms),collapse = " + "))
    form <- as.formula(form)
  } else {
    lLevels <- paste0("level",levels(data$level))#hierFormula$Get("name",filterFun=function(x) x$level==2)
    formulaTerms <- formTerms <- recharge <- list()
    factorTerms <- names(data)[which(unlist(lapply(data,function(x) inherits(x,"factor"))))]
    for(j in lLevels){
      newForm <- newFormulas(hierFormula[[j]]$formula, length(hierStates$Get("state",filterFun=data.tree::isLeaf)))
      #if(j==lLevels[1]) if(!attr(stats::terms(hierFormula[[j]]$formula),"intercept")) stop("hierFormula$",j,"$formula must include an intercept term")
      formulaTerms[[j]] <- newForm$formterms
      recharge[[j]] <- newForm$recharge
      formTerms[[j]] <- NULL
      if(length(formulaTerms[[j]])){
        if(any(grepl("level",formulaTerms[[j]]))) stop("hierFormula$",j,"$formula cannot include 'level'")
        if(any(formulaTerms[[j]] %in% factorTerms)){
          mm<-model.matrix(as.formula(paste0("~",as.character(hierFormula[[j]]$formula)[2])),data)
          as <- attr(mm,"assign")
          colmm <- colnames(mm)
          for(jj in 1:length(colmm)){
            cterms <- unlist(strsplit(formulaTerms[[j]][as[jj]],":"))
            lInd <- unlist(strsplit(colmm[jj],":"))[which(cterms %in% factorTerms)]
            for(kk in cterms[which(cterms %in% factorTerms)]){
              lInd <- gsub(kk,"",lInd)
            }
            names(lInd) <- cterms[which(cterms %in% factorTerms)]
            for(kk in cterms[which(cterms %in% factorTerms)]){
              colmm[jj] <- gsub(paste0(kk,lInd[kk]),paste0("(",kk,"=='",lInd[kk],"')"),colmm[jj])
            }
            colmm[jj] <- paste0("I((level=='",gsub("level","",j),"')*",colmm[jj],")")
            colmm[jj] <- gsub(":","*",colmm[jj])
            colmm[jj] <- gsub("*(Intercept)","*1",colmm[jj],fixed=TRUE)
          }
          formTerms[[j]] <- c(formTerms[[j]],colmm)
        } else {
          for(k in 1:length(formulaTerms[[j]])){
            formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*",formulaTerms[[j]][k],")"))
          }
        }
      }
      if(!is.null(recharge[[j]])){
        #formulaTerms[[j]] <- c(formulaTerms[[j]],paste0("recharge",gsub("level","",j)))
        #formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*","recharge",gsub("level","",j),")"))
        #data[[paste0("recharge",gsub("level","",j))]] <- rep(0,nrow(data))
        formulaTerms[[j]] <- c(formulaTerms[[j]],"recharge")
        formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*","recharge)"))
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
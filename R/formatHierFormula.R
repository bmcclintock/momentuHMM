#' @importFrom stats terms
#' @importFrom data.tree Get
formatHierFormula <- function(data,hierFormula){
  if(is.null(data)){
    lLevels <- hierFormula$Get("name",filterFun=function(x) x$level==2)
    formulaTerms <- list()
    for(j in lLevels){
      #if(j==lLevels[1]) if(!attr(stats::terms(hierFormula[[j]]$formula),"intercept")) stop("hierFormula$",j,"$formula must include an intercept term")
      formulaTerms[[j]] <- attr(stats::terms(hierFormula[[j]]$formula),"term.labels")
      if(length(formulaTerms[[j]])){
        if(any(grepl("level",formulaTerms[[j]]))) stop("hierFormula$",j,"$formula cannot include 'level'")
        formulaTerms[[j]] <- paste0("I((level=='",gsub("level","",j),"')*",formulaTerms[[j]],")")
      }
    }
    form <- "~ 0 + level"
    if(length(unlist(formulaTerms))) form <- paste0(form," + ",paste0(unlist(formulaTerms),collapse = " + "))
    form <- as.formula(form)
  } else {
    lLevels <- paste0("level",levels(data$level))#hierFormula$Get("name",filterFun=function(x) x$level==2)
    formulaTerms <- formTerms <- list()
    factorTerms <- names(data)[which(unlist(lapply(data,function(x) inherits(x,"factor"))))]
    for(j in lLevels){
      #if(j==lLevels[1]) if(!attr(stats::terms(hierFormula[[j]]$formula),"intercept")) stop("hierFormula$",j,"$formula must include an intercept term")
      formulaTerms[[j]] <- attr(stats::terms(hierFormula[[j]]$formula),"term.labels")
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
  return(form)
}
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
        if("level" %in% formulaTerms[[j]]) stop("hierFormula$",j,"$formula cannot include 'level'")
        formulaTerms[[j]] <- paste0("I((level=='",gsub("level","",j),"')*",formulaTerms[[j]],")")
      }
    }
    form <- "~ 0 + level"
    if(length(unlist(formulaTerms))) form <- paste0(form," + ",paste0(unlist(formulaTerms),collapse = " + "))
    form <- as.formula(form)
  } else {
    lLevels <- hierFormula$Get("name",filterFun=function(x) x$level==2)
    formulaTerms <- formTerms <- list()
    factorTerms <- names(data)[which(unlist(lapply(data,function(x) inherits(x,"factor"))))]
    for(j in lLevels){
      #if(j==lLevels[1]) if(!attr(stats::terms(hierFormula[[j]]$formula),"intercept")) stop("hierFormula$",j,"$formula must include an intercept term")
      formulaTerms[[j]] <- attr(stats::terms(hierFormula[[j]]$formula),"term.labels")
      formTerms[[j]] <- NULL
      if(length(formulaTerms[[j]])){
        if("level" %in% formulaTerms[[j]]) stop("hierFormula$",j,"$formula cannot include 'level'")
        for(k in 1:length(formulaTerms[[j]])){
          if(formulaTerms[[j]][k] %in% factorTerms) formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*(",formulaTerms[[j]][k],"=='",levels(data[[formulaTerms[[j]][k]]]),"'))"))
          else formTerms[[j]] <- c(formTerms[[j]],paste0("I((level=='",gsub("level","",j),"')*",formulaTerms[[j]][k],")"))
        }
      }
    }
    factorLevels <- which(unlist(lapply(formulaTerms,function(x) any(x %in% factorTerms))))
    if(any(factorLevels)){
      if(length(factorLevels)==length(lLevels)) form <- "~ 0"
      else form <- paste0("~ 0 + ",paste0("I((level=='",gsub("level","",lLevels[-factorLevels]),"')*1)",collapse=" + "))
    } else form <- "~ 0 + level"
    if(length(unlist(formulaTerms))) form <- paste0(form," + ",paste0(unlist(formTerms),collapse = " + "))
    form <- as.formula(form)
  }
  return(form)
}
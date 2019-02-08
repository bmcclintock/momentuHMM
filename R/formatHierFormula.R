#' @importFrom stats terms
#' @importFrom data.tree Get
formatHierFormula <- function(hierFormula){
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
  as.formula(form)
}
#' @importFrom survival untangle.specials
cosinorCos<-function(x,period){
  cos(2*pi*x/period)
}
cosinorSin<-function(x,period){
  sin(2*pi*x/period)
}

stateFormulas<-function(formula,nbStates,spec="state"){
  
  Terms <- terms(formula, specials = c(paste0(spec,1:nbStates),"cosinor"))
  
  stateFormula<-list()
  if(length(unlist(attr(Terms,"specials")))){
    varnames <- attr(Terms,"term.labels")
    mainpart <- varnames
    if(!is.null(unlist(attr(Terms, "specials")))){
      cosInd <- survival::untangle.specials(Terms,"cosinor",order=1:10)$terms
      stateInd <- numeric()
      for(j in 1:nbStates){
        tmpInd <- survival::untangle.specials(Terms,paste0(spec,j),order=1)$terms
        if(length(tmpInd)) stateInd<-c(stateInd,tmpInd)
      }
      if(length(cosInd) | length(stateInd)){
        mainpart <- varnames[-c(cosInd,stateInd)]
      }
      for(j in varnames[cosInd])
        mainpart<-c(mainpart,paste0(gsub("cosinor","cosinorCos",j)),paste0(gsub("cosinor","cosinorSin",j)))
    }
    
  
    for(j in 1:nbStates){
      tmplabs<-attr(Terms,"term.labels")[attr(Terms,"specials")[[paste0(spec,j)]]]
      if(length(tmplabs)){
        tmp<- terms(as.formula(paste("~",substr(tmplabs,nchar(paste0(spec,j))+1,nchar(tmplabs)))),specials="cosinor")
      
        tmpnames<-attr(tmp,"term.labels")
        mp<-tmpnames
        if(!is.null(unlist(attr(tmp,"specials")))){
          cosInd <- survival::untangle.specials(tmp,"cosinor",order=1:10)$terms
          mp <- c(tmpnames[-cosInd])
          for(i in tmpnames[cosInd])
            mp<-c(mp,paste0(gsub("cosinor","cosinorCos",i)),paste0(gsub("cosinor","cosinorSin",i)))
        }
      } else {
        tmp <- Terms
        mp <- character()
      }
      stateFormula[[j]]<-as.formula(paste("~",paste(c(attr(tmp,"intercept"),mainpart,mp),collapse = " + "),collapse=" + "))
    }
  } else {
    for(j in 1:nbStates){
      stateFormula[[j]] <- formula
    }
  }
  stateFormula
}
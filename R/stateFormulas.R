#' @importFrom survival untangle.specials
#' @importFrom prodlim strip.terms
cosinorCos<-function(x,period){
  cos(2*pi*x/period)
}
cosinorSin<-function(x,period){
  sin(2*pi*x/period)
}

#angleStrength<-function(angle,strength){
#  tt <- terms(f,specials=c("angleStrength"))
#  st <- strip.terms(tt,specials=c("angleStrength"),arguments=list(angleStrength=list("strength"=1)))
#}

stateFormulas<-function(formula,nbStates,spec="state",angleMean=FALSE){
  
  Terms <- terms(formula, specials = c(paste0(spec,1:nbStates),"cosinor","angleStrength"))
  
  stateFormula<-list()
  if(length(unlist(attr(Terms,"specials"))) | angleMean){
    varnames <- attr(Terms,"term.labels")
    mainpart <- varnames
    cosInd <- survival::untangle.specials(Terms,"cosinor",order=1:10)$terms
    if(length(cosInd) & angleMean) stop("cosinor models are not supported for angle means")
    angInd <- survival::untangle.specials(Terms,"angleStrength",order=1:10)$terms
    if(length(angInd) & !angleMean) stop("angleStrength models are only allowed for angle means")
    stateInd <- numeric()
    for(j in 1:nbStates){
      tmpInd <- survival::untangle.specials(Terms,paste0(spec,j),order=1)$terms
      if(length(tmpInd)) stateInd<-c(stateInd,tmpInd)
    }
    if(length(cosInd) | length(angInd) | length(stateInd)){
      mainpart <- varnames[-c(cosInd,angInd,stateInd)]
    }
    if(angleMean & length(mainpart)){
      tmpmainpart <- mainpart
      mainpart <- character()
      for(j in 1:length(tmpmainpart)){
        mainpart <- c(mainpart,paste0(c("sin","cos"),"(",tmpmainpart[j],")"))
      }
    }
    for(j in varnames[cosInd])
      mainpart<-c(mainpart,paste0(gsub("cosinor","cosinorCos",j)),paste0(gsub("cosinor","cosinorSin",j)))
    if(length(angInd)){
      stmp <- prodlim::strip.terms(Terms[attr(Terms,"specials")$angleStrength],specials="angleStrength",arguments=list(angleStrength=list("strength"=1)))
      for(jj in attr(stmp,"term.labels")){
        mainpart<-c(mainpart,paste0(attr(stmp,"stripped.arguments")$angleStrength[[jj]]$strength,":",c("sin","cos"),"(",jj,")"))
      }
    }
    
    for(j in 1:nbStates){
      tmplabs<-attr(Terms,"term.labels")[attr(Terms,"specials")[[paste0(spec,j)]]]
      if(length(tmplabs)){
        tmp<- terms(as.formula(paste("~",substr(tmplabs,nchar(paste0(spec,j))+1,nchar(tmplabs)),collapse="+")),specials=c("cosinor","angleStrength"))
      
        tmpnames<-attr(tmp,"term.labels")
        mp<-tmpnames
        if(!is.null(unlist(attr(tmp,"specials"))) | angleMean){
          cosInd <- survival::untangle.specials(tmp,"cosinor",order=1:10)$terms
          if(length(cosInd) & angleMean) stop("cosinor models are not supported for angle means")
          angInd <- survival::untangle.specials(tmp,"angleStrength",order=1:10)$terms
          if(length(angInd) & !angleMean) stop("angleStrength models are only allowed for angle means")
          if(length(cosInd) | length(angInd)){
            mp <- c(tmpnames[-c(cosInd,angInd)])
          }
          if(angleMean & length(mp)){
            tmpmp <- mp
            mp <- character()
            for(jj in 1:length(tmpmp)){
              mp <- c(mp,paste0(c("sin","cos"),"(",tmpmp[jj],")"))
            }
          }
          for(i in tmpnames[cosInd])
            mp<-c(mp,paste0(gsub("cosinor","cosinorCos",i)),paste0(gsub("cosinor","cosinorSin",i)))
          if(length(angInd)){
            stmp <- prodlim::strip.terms(tmp[attr(tmp,"specials")$angleStrength],specials="angleStrength",arguments=list(angleStrength=list("strength"=1)))
            for(jj in attr(stmp,"term.labels")){
              mp<-c(mp,paste0(attr(stmp,"stripped.arguments")$angleStrength[[jj]]$strength,":",c("sin","cos"),"(",jj,")"))
            }
          }
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
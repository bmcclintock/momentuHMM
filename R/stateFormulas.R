#' @importFrom survival untangle.specials
#' @importFrom prodlim strip.terms
cosinorCos<-function(x,period){
  cos(2*pi*x/period)
}
cosinorSin<-function(x,period){
  sin(2*pi*x/period)
}

#angleFormula<-function(angle,strength){
#  tt <- terms(f,specials=c("angleFormula"))
#  st <- strip.terms(tt,specials=c("angleFormula"),arguments=list(angleFormula=list("strength"=1)))
#}

recharge<-function(g0=~1, theta){
  if(!is.formula(g0)) stop("recharge function must include a formula for g0")
  if(!is.formula(theta)) stop("recharge function must include a formula for theta")
  return(list(g0=g0,theta=theta))
}

stateFormulas<-function(formula,nbStates,spec="state",angleMean=FALSE,data=NULL){
  
  Terms <- terms(formula, specials = c(paste0(spec,1:nbStates),"cosinor","angleFormula","recharge"))
  if(any(grepl("angleStrength\\(",attr(Terms,"term.labels")))) stop("'angleStrength' is defunct in momentuHMM >=1.4.2. Please use 'angleFormula' instead")
  if(any(attr(Terms,"order")>1)){
    if(any(grepl("angleFormula\\(",attr(Terms,"term.labels")[attr(Terms,"order")>1]))) stop("interactions with angleFormula are not allowed")
  }
  
  stateFormula<-list()
  if(length(unlist(attr(Terms,"specials"))) | angleMean){
    varnames <- attr(Terms,"term.labels")
    mainpart <- varnames
    cosInd <- survival::untangle.specials(Terms,"cosinor",order=1:10)$terms
    if(length(cosInd) & angleMean) stop("cosinor models are not supported for angle means")
    angInd <- survival::untangle.specials(Terms,"angleFormula",order=1:10)$terms
    if(length(angInd) & !angleMean) stop("angleFormula models are only allowed for angle means")
    recInd <- survival::untangle.specials(Terms,"recharge",order=1:10)$terms
    if(length(recInd)>1) stop("only a single recharge model is permitted")
    stateInd <- numeric()
    for(j in 1:nbStates){
      tmpInd <- survival::untangle.specials(Terms,paste0(spec,j),order=1)$terms
      if(length(tmpInd)) stateInd<-c(stateInd,tmpInd)
    }
    if(length(cosInd) | length(angInd) | length(recInd) | length(stateInd)){
      mainpart <- varnames[-c(cosInd,angInd,recInd,stateInd)]
    }
    if(angleMean & length(mainpart)){
      tmpmainpart <- mainpart
      mainpart <- character()
      if(any(grepl("cos",tmpmainpart)) | any(grepl("sin",tmpmainpart))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
      for(j in 1:length(tmpmainpart)){
        mainpart <- c(mainpart,paste0(c("sin","cos"),"(",tmpmainpart[j],")"))
      }
    }
    for(j in varnames[cosInd])
      mainpart<-c(mainpart,paste0(gsub("cosinor","cosinorCos",j)),paste0(gsub("cosinor","cosinorSin",j)))
    for(j in varnames[recInd])
      mainpart<-c(mainpart,paste0(gsub(")\\s*$","",gsub("recharge\\(","",j))))
    if(length(angInd)){
      stmp <- prodlim::strip.terms(Terms[attr(Terms,"specials")$angleFormula],specials="angleFormula",arguments=list(angleFormula=list("strength"=NULL,"by"=NULL)))
      if(any(grepl("cos",attr(stmp,"term.labels"))) | any(grepl("sin",attr(stmp,"term.labels")))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
      for(jj in attr(stmp,"term.labels")){
        if(is.null(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength)){
          tmpForm <- ~ - 1
          strengthInd <- FALSE
        }
        else {
          if(!is.na(suppressWarnings(as.numeric((attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))))) stop("angleFormula has invalid strength argument")
          tmpForm <- as.formula(paste0("~-1+",attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))
          if(any(attr(terms(tmpForm),"order")>1)) stop("angleFormula strength argument for ",jj," cannot include term interactions; use the 'by' argument")
          strengthInd <- TRUE
        }
        group <- attr(stmp,"stripped.arguments")$angleFormula[[jj]]$by
        if(!is.null(data)){
          DMterms <- attr(terms(tmpForm),"term.labels")
          factorterms<-names(data)[unlist(lapply(data,is.factor))]
          factorcovs<-paste0(rep(factorterms,times=unlist(lapply(data[factorterms],nlevels))),unlist(lapply(data[factorterms],levels)))
          for(cov in DMterms){
            form<-formula(paste("~",cov))
            varform<-all.vars(form)
            if(any(varform %in% factorcovs)){
              factorvar<-factorcovs %in% varform
              tmpcov<-rep(factorterms,times=unlist(lapply(data[factorterms],nlevels)))[which(factorvar)]
              tmpcovj <- cov
              for(j in 1:length(tmpcov)){
                tmpcovj <- gsub(factorcovs[factorvar][j],tmpcov[j],tmpcovj)
              }
              form <- formula(paste("~ 0 + ",tmpcovj))
            }
            if(strengthInd){
              if(any(model.matrix(form,data)<0)) stop(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength," must be >=0 in order to be used in angleFormula")
              if(any(unlist(lapply(data[all.vars(form)],function(x) inherits(x,"factor"))))) stop("angleFormula strength argument cannot be a factor; use the 'by' argument for factors")
            }
          }
          if(!is.null(group)){
            form<-formula(paste("~",group))
            varform<-all.vars(form)
            if(any(varform %in% factorcovs)){
              factorvar<-factorcovs %in% varform
              varform<-rep(factorterms,times=unlist(lapply(data[factorterms],nlevels)))[which(factorvar)]
            }
            if(any(!(varform %in% names(data)))) stop("angleFormula 'by' argument ",varform[which(!(varform %in% names(data)))]," not found in data")
            if(any(!unlist(lapply(data[varform],function(x) inherits(x,"factor"))))) stop("angleFormula 'by' argument must be of class factor")
          }
        }
        
        if(!is.null(group)) group <- paste0(group,":")
        mainpart<-c(mainpart,paste0(group,ifelse(strengthInd,paste0(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength,":"),""),c("sin","cos"),"(",jj,")"))
      }
    }
    
    for(j in 1:nbStates){
      tmplabs<-attr(Terms,"term.labels")[attr(Terms,"specials")[[paste0(spec,j)]]]
      if(length(tmplabs)){
        tmp<- terms(as.formula(paste("~",substr(tmplabs,nchar(paste0(spec,j))+1,nchar(tmplabs)),collapse="+")),specials=c("cosinor","angleFormula","recharge"))
        
        tmpnames<-attr(tmp,"term.labels")
        if(any(grepl("angleStrength\\(",tmpnames))) stop("'angleStrength' is defunct in momentuHMM >=1.4.2. Please use 'angleFormula' instead")
        if(any(attr(tmp,"order")>1)){
          if(any(grepl("angleFormula\\(",tmpnames[attr(tmp,"order")>1]))) stop("interactions with angleFormula are not allowed")
        }
        mp<-tmpnames
        if(!is.null(unlist(attr(tmp,"specials"))) | angleMean){
          cosInd <- survival::untangle.specials(tmp,"cosinor",order=1:10)$terms
          if(length(cosInd) & angleMean) stop("cosinor models are not supported for angle means")
          angInd <- survival::untangle.specials(tmp,"angleFormula",order=1:10)$terms
          if(length(angInd) & !angleMean) stop("angleFormula models are only allowed for angle means")
          recInd <- survival::untangle.specials(tmp,"recharge",order=1:10)$terms
          if(length(recInd)) stop("recharge models cannot be state-dependent")
          if(length(cosInd) | length(angInd)){
            mp <- c(tmpnames[-c(cosInd,angInd)])
          }
          if(angleMean & length(mp)){
            tmpmp <- mp
            mp <- character()
            if(any(grepl("cos",tmpmp)) | any(grepl("sin",tmpmp))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
            for(jj in 1:length(tmpmp)){
              mp <- c(mp,paste0(c("sin","cos"),"(",tmpmp[jj],")"))
            }
          }
          for(i in tmpnames[cosInd])
            mp<-c(mp,paste0(gsub("cosinor","cosinorCos",i)),paste0(gsub("cosinor","cosinorSin",i)))
          for(i in tmpnames[recInd])
            mp<-c(mp,paste0(gsub(")\\s*$","",gsub("recharge\\(","",i))))
          if(length(angInd)){
            stmp <- prodlim::strip.terms(tmp[attr(tmp,"specials")$angleFormula],specials="angleFormula",arguments=list(angleFormula=list("strength"=NULL,"by"=NULL)))
            if(any(grepl("cos",attr(stmp,"term.labels"))) | any(grepl("sin",attr(stmp,"term.labels")))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
            for(jj in attr(stmp,"term.labels")){
              if(is.null(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength)){
                tmpForm <- ~ - 1
                strengthInd <- FALSE
              }
              else {
                if(!is.na(suppressWarnings(as.numeric((attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))))) stop("angleFormula has invalid strength argument")
                tmpForm <- as.formula(paste0("~-1+",attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))
                if(any(attr(terms(tmpForm),"order")>1)) stop("angleFormula strength argument cannot include term interactions; use the 'by' argument")
                strengthInd <- TRUE
              }
              group <- attr(stmp,"stripped.arguments")$angleFormula[[jj]]$by
              if(!is.null(data)){
                if(strengthInd){
                  if(any(model.matrix(tmpForm,data)<0)) stop(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength," must be >=0 in order to be used in angleFormula")
                  if(inherits(data[[attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength]],"factor")) stop("angleFormula strength argument cannot be a factor; use the 'by' argument for factors")
                }
                if(!is.null(group)){
                  varform <- all.vars(as.formula(paste0("~",group)))
                  if(any(!(varform %in% names(data)))) stop("angleFormula 'by' argument ",varform[which(!(varform %in% names(data)))]," not found in data")
                  if(any(!unlist(lapply(data[varform],function(x) inherits(x,"factor"))))) stop("angleFormula 'by' argument must be of class factor")
                }
              }
              if(!is.null(group)) group <- paste0(group,":")
              mp<-c(mp,paste0(group,ifelse(strengthInd,paste0(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength,":"),""),c("sin","cos"),"(",jj,")"))
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

newFormulas<-function(formula,nbStates,hierarchical=FALSE)
{
  stateForms<- terms(formula, specials = c(paste0(rep(c("state","toState"),each=nbStates),1:nbStates),"recharge"))
  newformula<-formula
  formulaStates <- vector('list',nbStates*(nbStates-1))
  formulaStates[1:(nbStates*(nbStates-1))] <- list(newformula)
  formterms<-attr(terms.formula(newformula),"term.labels")
  recharge <- NULL
  
  if(nbStates>1){
    if(length(unlist(attr(stateForms,"specials")))){
      #newForm<-attr(stateForms,"term.labels")[-unlist(attr(stateForms,"specials"))]
      newForm <- attr(stateForms,"term.labels")[-which(colSums(attr(stateForms,"factors")[unlist(attr(stateForms,"specials")),,drop=FALSE])>0)]#attr(drop.terms(stateForms,which(attr(stateForms,"factors")[unlist(attr(stateForms,"specials")),]>0)),"term.labels")
      for(i in 1:nbStates){
        if(!is.null(attr(stateForms,"specials")[[paste0("state",i)]])){
          for(j in 1:(nbStates-1)){
            #newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("state",i)]]]))
            stateFact <- attr(stateForms,"factors")[attr(stateForms,"specials")[[paste0("state",i)]],,drop=FALSE]
            newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),colnames(stateFact)[which(stateFact>0)]))
          }
        }
        if(!is.null(attr(stateForms,"specials")[[paste0("toState",i)]])){
          betaInd<-matrix(0,nbStates,nbStates,byrow=TRUE)
          diag(betaInd)<-NA
          betaInd[!is.na(betaInd)] <- seq(1:(nbStates*(nbStates-1)))
          betaInd<-t(betaInd)[,i]
          betaInd<-betaInd[!is.na(betaInd)]
          for(j in betaInd){
            #newForm<-c(newForm,gsub(paste0("toState",i),paste0("betaCol",j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("toState",i)]]]))
            stateFact <- attr(stateForms,"factors")[attr(stateForms,"specials")[[paste0("toState",i)]],,drop=FALSE]
            newForm<-c(newForm,gsub(paste0("toState",i),paste0("betaCol",j),colnames(stateFact)[which(stateFact>0)]))
          }
        }
      }
      if(!is.null(attr(stateForms,"specials")$recharge)){
        reSub <- rownames(attr(stateForms,"factors"))[attr(stateForms,"specials")$recharge]
        recharge <- lapply(reSub,function(x) eval(parse(text=x)))
        #recharge <- stateFormulas(as.formula(paste("~",paste(attr(stateForms,"term.labels")[which(grepl("recharge",attr(stateForms,"term.labels")))]))),1)[[1]]
        hierInd <- FALSE
        for(k in 1:length(reSub)){
          #if(!hierarchical){
            reFact <- attr(stateForms,"factors")[attr(stateForms,"specials")$recharge[k],,drop=FALSE]
            iLevel <- gsub("\") * 1)","",gsub("I((level == \"","",names(which(attr(stateForms,"factors")[,colnames(reFact)[which(reFact>0)]]==2)),fixed=TRUE),fixed=TRUE)
            if(length(iLevel)){
              hierInd <- TRUE
              if(!hierarchical){
                newForm <- c(newForm,gsub(reSub[k],paste0("recharge",iLevel),colnames(reFact)[which(reFact>0)],fixed=TRUE))
              }
              names(recharge)[k] <- paste0("level",iLevel)
            } else {
              if(!hierarchical){
                newForm <- c(newForm,gsub(reSub[k],"recharge",colnames(reFact)[which(reFact>0)],fixed=TRUE))
              }
            }
          #}
        }
        if(length(newForm)){
          newformula<-as.formula(paste("~",ifelse(attr(terms(formula),"intercept"),"","0+"),paste(newForm,collapse="+")))
        } else {
          newformula <- ~1
        }
        if(!hierInd) recharge <- recharge[[1]]
      } else {
        newformula<-as.formula(paste("~",paste(newForm,collapse="+")))
      }
    }
    formulaStates<-stateFormulas(newformula,nbStates*(nbStates-1),spec="betaCol")
    if(length(unlist(attr(terms(newformula, specials = c(paste0("betaCol",1:(nbStates*(nbStates-1))),"cosinor")),"specials")))){
      allTerms<-unlist(lapply(formulaStates,function(x) attr(terms(x),"term.labels")))
      newformula<-as.formula(paste("~",paste(allTerms,collapse="+")))
      formterms<-attr(terms.formula(newformula),"term.labels")
    } else if(is.null(recharge)) {
      formterms<-attr(terms.formula(newformula),"term.labels")
      newformula<-formula
    } else {
      formterms<-attr(terms.formula(newformula),"term.labels")
    }
    if(hierarchical){
      recharge <- expandRechargeFormulas(recharge)
    }
  }  
  return(list(formulaStates=formulaStates,formterms=formterms,newformula=newformula,recharge=recharge))
}
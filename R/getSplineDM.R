#' @importFrom stats predict
# #' @importFrom qdapRegex rm_between
getSplineDM<-function(distnames,DM,m,covs){
  splineInd<-list()
  splineCovs<-list()
  newDM<-DM
  newcovs<-covs
  nbStates<-length(m$stateNames)
  for(i in distnames){
    splineInd[[i]]<-list()
    splineCovs[[i]]<-list()
    if(!is.null(DM[[i]])){
      if(is.list(DM[[i]])){
        newDM[[i]]<-list()
        for(j in names(DM[[i]])){
          splineInd[[i]][[j]]<-FALSE
          splineCovs[[i]][[j]]<-character()
          formulaStates<-stateFormulas(DM[[i]][[j]],nbStates,angleMean=(j=="mean" & !isFALSE(m$conditions$circularAngleMean[[i]])))
          tmpDM<-list()
          for(state in 1:nbStates){
            tmpDM[[state]]<-character()
            Terms<-stats::terms(formulaStates[[state]],specials=splineList)
            factors<-attr(Terms,"factors")
            specials<-rownames(factors)[unlist(attr(Terms,"specials"))]
            for(k in rownames(factors)){
              if(k %in% specials){
                splineInd[[i]][[j]]<-TRUE
                splineCovs[[i]][[j]]<-unique(c(splineCovs[[i]][[j]],all.vars(stats::as.formula(paste0("~",k)))))
                if (!requireNamespace("qdapRegex", quietly = TRUE)) {
                  stop("Package \"qdapRegex\" needed for this function to work. Please install it.",
                       call. = FALSE)
                }
                splineExpr<-qdapRegex::rm_between(k, "(", ",", extract=TRUE)[[1]]
                sp<-eval(substitute(eval(parse(text=k))),m$data,parent.frame())
                tmpcovs<-predict(sp,eval(substitute(eval(parse(text=splineExpr))),covs,parent.frame()))
                tmp<-colnames(stats::model.matrix(stats::as.formula(paste0("~",k)),m$data)[,-1])
                tmp<-gsub("[()]","",tmp)
                tmp<-gsub(" ","_",tmp)
                tmp<-gsub(",","_",tmp)
                tmp<-gsub("=","_",tmp)
                tmp<-gsub("\\*","_",tmp)
                tmp<-gsub("\\^","_",tmp)
                colnames(tmpcovs)<-tmp
                newcovs<-cbind(newcovs,tmpcovs)
                for(l in colnames(factors)){
                  if(factors[k,l]){
                    newcols<-stats::model.matrix(stats::as.formula(paste0("~",l)),m$data)[,-1]
                    tmp<-colnames(newcols)
                    tmp<-gsub("[()]","",tmp)
                    tmp<-gsub(" ","_",tmp)
                    tmp<-gsub(",","_",tmp)
                    tmp<-gsub("=","_",tmp)
                    tmp<-gsub("\\*","_",tmp)
                    tmp<-gsub("\\^","_",tmp)
                    tmpDM[[state]]<-c(tmpDM[[state]],tmp)
                  }
                }
              } else {
                if(length(specials)) {
                  tmpspec<-specials
                  tmpspec<-gsub("(","\\(",tmpspec,fixed=TRUE)
                  tmpspec<-gsub(")","\\)",tmpspec,fixed=TRUE)
                  lfact <- grep(paste(tmpspec,collapse="|"),colnames(factors), value=TRUE,invert=TRUE)
                } else lfact <- colnames(factors)
                for(l in lfact){
                  if(factors[k,l]){
                    tmpDM[[state]]<-c(tmpDM[[state]],l)
                  }
                }
              }
            }
          }
          if(!splineInd[[i]][[j]]) newDM[[i]][[j]] <- DM[[i]][[j]]
          else {
            tmpterms<-character()
            for(state in 1:nbStates){
              if(length(tmpDM[[state]]))
                tmpterms<-c(tmpterms,paste0("state",state,"(",paste0(tmpDM[[state]],collapse="+"),")"))
            }
            newDM[[i]][[j]] <- stats::as.formula(paste0("~",paste0(attr(DM[[i]][[j]],"intercept"),tmpterms,collapse="+")))
          }
        }
      }
    }
  }
  return(list(DM=newDM,covs=newcovs))
}

getSplineFormula<-function(formula,data,covs){
  newcovs<-covs
  splineInd<-FALSE
  splineCovs<-character()
  newformula<-character()
  Terms<-stats::terms(formula,specials=splineList)
  factors<-attr(Terms,"factors")
  specials<-rownames(factors)[unlist(attr(Terms,"specials"))]
  for(k in rownames(factors)){
    if(k %in% specials){
      splineInd<-TRUE
      splineCovs<-unique(c(splineCovs,all.vars(stats::as.formula(paste0("~",k)))))
      if (!requireNamespace("qdapRegex", quietly = TRUE)) {
        stop("Package \"qdapRegex\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
      splineExpr<-qdapRegex::rm_between(k, "(", ",", extract=TRUE)[[1]]
      sp<-eval(substitute(eval(parse(text=k))),data,parent.frame())
      tmpcovs<-predict(sp,eval(substitute(eval(parse(text=splineExpr))),covs,parent.frame()))
      tmp<-colnames(stats::model.matrix(stats::as.formula(paste0("~",k)),data)[,-1])
      tmp<-gsub("[()]","",tmp)
      tmp<-gsub(" ","_",tmp)
      tmp<-gsub(",","_",tmp)
      tmp<-gsub("=","_",tmp)
      colnames(tmpcovs)<-tmp
      newcovs<-cbind(newcovs,tmpcovs)
      for(l in colnames(factors)){
        if(factors[k,l]){
          tmp<-colnames(stats::model.matrix(stats::as.formula(paste0("~",l)),data))[-1]
          tmp<-gsub("[()]","",tmp)
          tmp<-gsub(" ","_",tmp)
          tmp<-gsub(",","_",tmp)
          tmp<-gsub("=","_",tmp)
          newformula<-c(newformula,tmp)
        }
      }
    } else {
      if(length(specials)) {
        tmpspec<-specials
        tmpspec<-gsub("(","\\(",tmpspec,fixed=TRUE)
        tmpspec<-gsub(")","\\)",tmpspec,fixed=TRUE)
        lfact <- grep(paste(tmpspec,collapse="|"),colnames(factors), value=TRUE,invert=TRUE)
      } else lfact <- colnames(factors)
      for(l in lfact){
        if(factors[k,l]){
          newformula<-c(newformula,l)
        }
      }
    }
  }
  if(!splineInd) newformula <- formula
  else newformula <- stats::as.formula(paste0("~",paste0(unique(newformula),collapse="+")))

  return(list(formula=newformula,covs=newcovs))
}

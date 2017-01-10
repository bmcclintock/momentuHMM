getboundInd<-function(DM=NULL){
  if(is.null(DM)){
    Ind<-NULL
  } else {
    Ind<-apply(apply(DM,1,function(x) apply(unique(DM),1,function(y) all(y==x))),2,function(x) which(x))
  }
  Ind
}

n2wDM<-function(bounds,DM,par,cons,logitcons){
  Par<-par
  tmpbounds<-bounds
  if(!is.numeric(bounds)){
    tmpind<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)!=bounds
    bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(bounds)))
    par <- par - bounds[,1]*tmpind[,1]
  } else {
    tmpind<-matrix(FALSE,nrow(bounds),ncol(bounds))    
  }
  a<-bounds[,1]
  for(j in which(tmpind[,1])){
    for(i in which(!tmpind[,1])){
      if(grepl(paste0("[",i,"]"),tmpbounds[j,1],fixed=TRUE)){
        a[j]<-bounds[i,1]
      }
    }
  }
  b<-bounds[,2]
  for(j in which(tmpind[,2])){
    for(i in which(!tmpind[,2])){
      if(grepl(paste0("[",i,"]"),tmpbounds[j,2],fixed=TRUE)){
        b[j]<-bounds[i,2]
      }
    }
  }
  ind1<-which(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind2<-which(!abs(a- -pi)<1.e-6 & !abs(b - pi)<1.e-6)
  p<-numeric(length(ind1)+length(ind2))
  if(length(ind1)) p[ind1] <- (solve(unique(DM)[ind1,ind1],tan(par[ind1]/2))-logitcons[ind1])^(1/cons[ind1])
  par<-par[ind2]
  a<-a[ind2]
  b<-b[ind2]
  cons<-cons[ind2]
  logitcons<-logitcons[ind2]
  if(all(is.finite(a)) & all(is.infinite(b))){
    p[ind2]<-solve(unique(DM)[ind2,ind2],log(par-a))^(1/cons)
  } else if(all(is.finite(a)) & all(is.finite(b))){
    p[ind2]<-(solve(unique(DM)[ind2,ind2],logit((par-a)/(b-a)))-logitcons)^(1/cons)
  } else {
    ind21<-which(is.infinite(b))
    ind22<-which(is.finite(b))
    p<-numeric(length(ind1)+length(ind2))
    p[ind2][ind21]<-solve(unique(DM)[ind21,ind21],log(par[ind21]-a[ind21]))^(1/cons[ind21])
    p[ind2][ind22]<-(solve(unique(DM)[ind22,ind22],logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))))^(1/cons[ind22])
  }
  p
}     

w2nDM<-function(wpar,bounds,DM,DMind,cons,logitcons,nbObs,k=0){
  Par<-numeric(length(wpar))
  #tmpind<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)!=bounds
  tmpbounds<-bounds
  if(!is.numeric(bounds)){
    tmpind<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)!=bounds
    bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(bounds)))
  } else {
    tmpind<-matrix(FALSE,nrow(bounds),ncol(bounds))
  }
  Bndind<-which(tmpind)
  #bndind<-which(tmpind[boundInd])
  a<-bounds[,1]
  for(j in which(tmpind[,1])){
    for(i in which(!tmpind[,1])){
      if(grepl(paste0("[",i,"]"),tmpbounds[j,1],fixed=TRUE)){
        a[j]<-bounds[i,1]
      }
    }
  }
  b<-bounds[,2]
  for(j in which(tmpind[,2])){
    for(i in which(!tmpind[,2])){
      if(grepl(paste0("[",i,"]"),tmpbounds[j,2],fixed=TRUE)){
        b[j]<-bounds[i,2]
      }
    }
  }
  a_2<-a#a[boundInd]
  b_2<-b#b[boundInd]
  #ind1<-which(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  #ind2<-which(!abs(a- -pi)<1.e-6 & !abs(b - pi)<1.e-6)
  ind11<-which(abs(a_2- -pi)<1.e-6 & abs(b_2 - pi)<1.e-6)
  ind21<-which(!abs(a_2- -pi)<1.e-6 & !abs(b_2 - pi)<1.e-6)
  
  if(DMind) XB <- matrix((wpar^cons+logitcons)%*%t(DM),nrow(DM),1)
  else XB <- getXB(DM,nbObs,wpar,cons,logitcons)
  
  p<-matrix(0,length(ind11)+length(ind21),nbObs)
  #p[ind11] <- 2*atan((wpar[ind1]^cons[ind1]+logitcons[ind1])%*%t(DM[ind1,ind1]))
  
  if(length(ind11))
    p[ind11,] <- (2*atan(XB))[ind11,]
  
  #Par<-numeric(length(ind1)+length(ind2))
  #Par[ind1] <- 2*atan((wpar[ind1]^cons[ind1]+logitcons[ind1])%*%t(unique(DM)[ind1,ind1]))
    
  #wpar<-wpar[ind2]
  #cons<-cons[ind2]
  #logitcons<-logitcons[ind2]
  
  
  #if(all(is.finite(a[ind2])) & all(is.infinite(b[ind2]))){
  if(all(is.finite(a[ind21])) & all(is.infinite(b[ind21]))){
    #Par[ind2]<-exp((wpar^cons)%*%t(unique(DM)[ind2,ind2]))+a[ind2]
    p[ind21,]<-(exp(XB)+a_2)[ind21,,drop=FALSE]
  #} else if(all(is.finite(a[ind2])) & all(is.finite(b[ind2]))){
  } else if(all(is.finite(a[ind21])) & all(is.finite(b[ind21]))){
    #Par[ind2]<-(b[ind2]-a[ind2])*(inv.logit((wpar^cons+logitcons)%*%t(unique(DM)[ind2,ind2])))+a[ind2]
    p[ind21,] <- ((b_2-a_2)*inv.logit(XB)+a_2)[ind21,,drop=FALSE]
  } else {
    #ind12<-which(is.infinite(b[ind2]))
    #ind22<-which(is.finite(b[ind2]))
    #Par[ind2][ind12]=exp(wpar[ind12]^cons[ind12]%*%t(unique(DM)[ind12,ind12]))+a[ind2][ind12]
    #Par[ind2][ind22]=(b[ind2][ind22]-a[ind2][ind22])*(inv.logit((wpar[ind22]^cons[ind22])%*%t(unique(DM)[ind22,ind22])))+a[ind2][ind22]
    
    ind13<-which(is.infinite(b_2[ind21]))
    ind23<-which(is.finite(b_2[ind21]))
    p[ind21[ind13],]<-(exp(XB)+a_2)[ind21[ind13],,drop=FALSE]
    p[ind21[ind23],]<-((b_2-a_2)*inv.logit(XB)+a_2)[ind21[ind23],,drop=FALSE]
  }
  #Par[Bndind]<-Par[Bndind]+matrix(sapply(tmpbounds,function(x) eval(parse(text=x))),ncol=2)[Bndind]
  #Par[which(Par>b)]<-b[which(Par>b)]
  if(length(Bndind)) p[Bndind,]<-p[Bndind,]+matrix(sapply(tmpbounds[Bndind],function(x) eval(parse(text=x))),ncol=2)
  if(any(p>b_2)) p[which(p>b_2)]<-b_2[which(p>b_2)%%length(b_2)]
  
  if(any(p<a_2 | p>b_2))
    stop("Scaling error.")
  
  if(k) {
    p <- p[k]
    return(p)
  } else {
    if(DMind) p<-matrix(p,nrow(p),nbObs)
    return(list(p=p,Par=Par))
  }
}   
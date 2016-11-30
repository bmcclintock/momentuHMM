getboundInd<-function(DM=NULL){
  if(is.null(DM)){
    Ind<-NULL
  } else {
    Ind<-apply(apply(DM,1,function(x) apply(unique(DM),1,function(y) all(y==x))),2,function(x) which(x))
  }
  Ind
}

n2wDM<-function(bounds,DM,par,cons,logitcons=rep(0,length(par))){
  Par<-par
  tmpind<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)!=bounds
  tmpbounds<-bounds
  if(!is.numeric(bounds)){
    bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(bounds)))
    par <- par - bounds[,1]*tmpind[,1]
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
  if(all(is.finite(a)) & all(is.infinite(b))){
    p<-solve(unique(DM),log(par-a))^(1/cons)
  } else if(all(is.finite(a)) & all(is.finite(b))){
    p<-(solve(unique(DM),logit((par-a)/(b-a)))-logitcons)^(1/cons)
  } else {
    ind1<-which(is.infinite(b))
    ind2<-which(is.finite(b))
    p=numeric(length(ind1)+length(ind2))
    p[ind1]=solve(unique(DM)[ind1,ind1],log(par[ind1]-a[ind1]))^(1/cons[ind1])
    p[ind2]=(solve(unique(DM)[ind2,ind2],logit((par[ind2]-a[ind2])/(b[ind2]-a[ind2]))))^(1/cons[ind2])
  }
  p
}     

w2nDM<-function(wpar,bounds,DM,cons,boundInd,logitcons=rep(0,length(wpar)),k=0){
  Par<-numeric(length(wpar))
  tmpind<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)!=bounds
  tmpbounds<-bounds
  if(!is.numeric(bounds)){
    bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(bounds)))
  }
  Bndind<-which(tmpind)
  bndind<-which(tmpind[boundInd])
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
  a_2<-a[boundInd]
  b_2<-b[boundInd]
  if(all(is.finite(a)) & all(is.infinite(b))){
    Par<-exp((wpar^cons)%*%t(unique(DM)))+a
    p<-exp((wpar^cons)%*%t(DM))+a_2
  } else if(all(is.finite(a)) & all(is.finite(b))){
    Par<-(b-a)*(inv.logit((wpar^cons+logitcons)%*%t(unique(DM))))+a
    p<-(b_2-a_2)*(inv.logit((wpar^cons+logitcons)%*%t(DM)))+a_2
  } else {
    ind1<-which(is.infinite(b))
    ind2<-which(is.finite(b))
    Par=numeric(length(ind1)+length(ind2))
    Par[ind1]=exp(wpar[ind1]^cons[ind1]%*%t(unique(DM)[ind1,ind1]))+a[ind1]
    Par[ind2]=(b[ind2]-a[ind2])*(inv.logit((wpar[ind2]^cons[ind2])%*%t(unique(DM)[ind2,ind2])))+a[ind2]
    
    ind11<-which(is.infinite(b))
    ind21<-which(is.finite(b))
    p=numeric(length(ind11)+length(ind21))
    p[ind11]=exp(wpar[ind1]^cons[ind1]%*%t(DM[ind11,ind1]))+a_2[ind11]
    p[ind21]=(b_2[ind21]-a_2[ind21])*(inv.logit((wpar[ind2]^cons[ind2])%*%t(DM[ind21,ind2])))+a_2[ind21]
  }
  Par[Bndind]<-Par[Bndind]+matrix(sapply(tmpbounds,function(x) eval(parse(text=x))),ncol=2)[Bndind]
  Par[which(Par>b)]<-b[which(Par>b)]
  p[bndind]<-p[bndind]+matrix(sapply(tmpbounds,function(x) eval(parse(text=x))),ncol=2)[boundInd,][bndind]
  p[which(p>b_2)]<-b_2[which(p>b_2)]
  if(k) {
    p <- p[k]
    return(p)
  } else {
    return(list(p=p,Par=Par))
  }
}   
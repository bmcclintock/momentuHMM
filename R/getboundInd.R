getboundInd<-function(DM=NULL){
  Ind<-1
  if(is.null(DM)){
    Ind<-NULL
  } else if(length(dim(DM))){
    Ind<-apply(apply(DM,1,function(x) apply(unique(DM),1,function(y) all(y==x))),2,function(x) which(x))
  }
  Ind
}
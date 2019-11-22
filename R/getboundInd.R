getboundInd<-function(DM=NULL){
  Ind<-1
  if(is.null(DM)){
    Ind<-NULL
  } else if(length(dim(DM))){
    uniDM <- matrix(unique(DM),ncol=dim(DM)[2])
    Ind<-apply(matrix(apply(DM,1,function(x) apply(uniDM,1,function(y) all(y==x))),dim(uniDM)[1],dim(DM)[1]),2,function(x) which(x))
  }
  Ind
}
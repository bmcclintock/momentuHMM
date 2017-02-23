momentuHMMdists<-c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm')
angledists<-c('wrpcauchy','vm')
stepdists<-c('gamma','weibull','exp','lnorm')
singleParmdists<-c('exp','pois')
nonnegativedists<-c('gamma','weibull','exp','lnorm','pois')
zeroInflationdists<-c('gamma','weibull','exp','lnorm','beta')

getboundInd<-function(DM=NULL){
  if(is.null(DM)){
    Ind<-NULL
  } else {
    Ind<-apply(apply(DM,1,function(x) apply(unique(DM),1,function(y) all(y==x))),2,function(x) which(x))
  }
  Ind
}
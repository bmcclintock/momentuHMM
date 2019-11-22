
#' Parameters definition
#'
#' @param dist Named list indicating the probability distributions of the data streams. 
#' @param nbStates Number of states of the HMM.
#' @param estAngleMean Named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy').
#' @param zeroInflation Named list of logicals indicating whether the probability distributions of the data streams should be zero-inflated.
#' @param oneInflation Named list of logicals indicating whether the probability distributions of the data streams are one-inflated.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of linear regression formulas or a matrix.  
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. For example, for a 2-state model using the wrapped Cauchy ('wrpcauchy') distribution for 
#' a data stream named 'angle' with \code{estAngleMean$angle=TRUE)}, \code{userBounds=list(angle=matrix(c(-pi,-pi,-1,-1,pi,pi,1,1),4,2))} 
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' 
#' @return A list of:
#' \item{parSize}{Named list indicating the number of natural parameters of the data stream probability distributions.}
#' \item{bounds}{Named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream.}
#' \item{parNames}{Names of parameters of the probability distribution for each data stream.}
#' \item{Bndind}{Named list indicating whether \code{DM} is NULL with default parameter bounds for each data stream.}
#' 
#' @examples 
#' \dontrun{
#' pD<-momentuHMM:::parDef(list(step="gamma",angle="wrpcauchy"),
#'     nbStates=2,list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),
#'     list(step=FALSE,angle=FALSE),NULL,NULL)
#' }

parDef <- function(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,userBounds=NULL)
{
  distnames<-names(dist)
  parSize <- parNames <- bounds <- vector('list',length(dist))
  names(parSize) <- names(parNames) <- names(bounds) <- distnames

  for(i in distnames){
    if(grepl("cat",dist[[i]])){
      dimCat <- as.integer(gsub("cat","",dist[[i]]))
      if(is.na(dimCat)) stop("categorical distributions must be specified using paste0('cat',k), where k is the number of categories (e.g. 'cat3', 'cat12', etc.)")
      if(dimCat<2) stop("categorical distribution must have at least 2 categories")
      dist[[i]] <- "cat"
    }
    switch(dist[[i]],
           "bern"={
             parSize[[i]] <- 1
             tmpbounds <- matrix(rep(c(0,1),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]]<-"prob"
           },
           "beta"={
             parSize[[i]] <- 2 + zeroInflation[[i]] + oneInflation[[i]]
             tmpbounds <- matrix(rep(c(0,Inf),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]]<-c("shape1","shape2")
           },
           "cat"={
             parSize[[i]] <- dimCat-1
             tmpbounds <- matrix(rep(c(0,1),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]]<-paste0("prob",1:parSize[[i]])
           },
           "exp"={
             parSize[[i]] <- 1 + zeroInflation[[i]]
             tmpbounds <- matrix(rep(c(0,Inf),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]] <- c("rate")
           },
           "gamma"={
             parSize[[i]] <- 2 + zeroInflation[[i]]
             tmpbounds <- matrix(rep(c(0,Inf),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]]<-c("mean","sd")
           },
           "logis"={
             parSize[[i]] <- 2
             tmpbounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                 ncol=2,byrow=TRUE)
             parNames[[i]] <- c("location","scale")
           },
           "lnorm"={
             parSize[[i]] <- 2 + zeroInflation[[i]]
             tmpbounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates*(1+zeroInflation[[i]]))),
                                  ncol=2,byrow=TRUE)
             parNames[[i]] <- c("location","scale")
           },
           "negbinom"={
             parSize[[i]] <- 2 
             tmpbounds <- matrix(rep(c(0,Inf),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]]<-c("mu","size")
           },
           "norm"={
             parSize[[i]] <- 2
             tmpbounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                  ncol=2,byrow=TRUE)
             parNames[[i]] <- c("mean","sd")
           },
           "mvnorm2"={
             parSize[[i]] <- 2 + 3
             tmpbounds <- matrix(c(rep(rep(c(-Inf,Inf),nbStates),2),rep(c(0,Inf),nbStates),rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                 ncol=2,byrow=TRUE)
             parNames[[i]] <- c("mean.x","mean.y","sigma.x","sigma.xy","sigma.y")
           },
           "mvnorm3"={
             parSize[[i]] <- 3 + 6
             tmpbounds <- matrix(c(rep(rep(c(-Inf,Inf),nbStates),3),rep(c(0,Inf),nbStates),
                                                                    rep(rep(c(-Inf,Inf),nbStates),2),
                                                                    rep(c(0,Inf),nbStates),
                                                                    rep(c(-Inf,Inf),nbStates),
                                                                    rep(c(0,Inf),nbStates)),ncol=2,byrow=TRUE)
             parNames[[i]] <- c("mean.x","mean.y","mean.z","sigma.x","sigma.xy","sigma.xz","sigma.y","sigma.yz","sigma.z")
           },
           "rw_norm"={
             parSize[[i]] <- 2
             tmpbounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                 ncol=2,byrow=TRUE)
             parNames[[i]] <- c("mean","sd")
           },
           "rw_mvnorm2"={
             parSize[[i]] <- 2 + 3
             tmpbounds <- matrix(c(rep(rep(c(-Inf,Inf),nbStates),2),rep(c(0,Inf),nbStates),rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                 ncol=2,byrow=TRUE)
             parNames[[i]] <- c("mean.x","mean.y","sigma.x","sigma.xy","sigma.y")
           },
           "rw_mvnorm3"={
             parSize[[i]] <- 3 + 6
             tmpbounds <- matrix(c(rep(rep(c(-Inf,Inf),nbStates),3),rep(c(0,Inf),nbStates),
                                   rep(rep(c(-Inf,Inf),nbStates),2),
                                   rep(c(0,Inf),nbStates),
                                   rep(c(-Inf,Inf),nbStates),
                                   rep(c(0,Inf),nbStates)),ncol=2,byrow=TRUE)
             parNames[[i]] <- c("mean.x","mean.y","mean.z","sigma.x","sigma.xy","sigma.xz","sigma.y","sigma.yz","sigma.z")
           },
           "pois"={
             parSize[[i]] <- 1
             tmpbounds <- matrix(rep(c(0,Inf),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]]<-"lambda"
           },
          "t"={
             parSize[[i]] <- 2
             tmpbounds <- matrix(c(rep(c(0,Inf),nbStates),rep(c(-Inf,Inf),nbStates)),
                                 ncol=2,byrow=TRUE)
             parNames[[i]] <- c("df","ncp")
           },
           "vm"={
             if(estAngleMean[[i]]) { # if the angle mean is estimated
               parSize[[i]] <- 2
               if(is.matrix(DM[[i]])){
                 dm <- DM[[i]]
                 meanind<-unique(unlist(apply(dm[1:nbStates,,drop=FALSE],1,function(x) which(x!=0))))
                 sdind<-unique(unlist(apply(dm[nbStates+1:nbStates,,drop=FALSE],1,function(x) which(x!=0))))
                 if(any(intersect(meanind,sdind))) stop("'DM' for ",i," cannot have parameters in common for mean and concentration")
               }
               tmpbounds <- rbind(matrix(rep(c(-pi,pi),nbStates),ncol=2,byrow=TRUE),matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE))
               parNames[[i]] <- c("mean","concentration") 
            }
             else {
               parSize[[i]] <- 1
               tmpbounds <- matrix(rep(c(0,Inf),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
               parNames[[i]] <- c("concentration")
             }
           },
           "vmConsensus"={
             parSize[[i]] <- 2
             if(is.matrix(DM[[i]])){
               dm <- DM[[i]]
               meanind<-unique(unlist(apply(dm[1:nbStates,,drop=FALSE],1,function(x) which(x!=0))))
               sdind<-unique(unlist(apply(dm[nbStates+1:nbStates,,drop=FALSE],1,function(x) which(x!=0))))
               if(any(intersect(meanind,sdind))) stop("'DM' for ",i," cannot have parameters in common for mean and kappa")
             }
             tmpbounds <- rbind(matrix(rep(c(-pi,pi),nbStates),ncol=2,byrow=TRUE),matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE))
             parNames[[i]] <- c("mean","kappa") 
           },
           "weibull"={
             parSize[[i]] <- 2 + zeroInflation[[i]]
             tmpbounds <- matrix(rep(c(0,Inf),parSize[[i]] * nbStates),ncol=2,byrow=TRUE)
             parNames[[i]]<-c("shape","scale")
           },
           "wrpcauchy"={
             if(estAngleMean[[i]]) {
               parSize[[i]] <- 2
               if(is.matrix(DM[[i]])){
                 dm <- DM[[i]]
                 meanind<-unique(unlist(apply(dm[1:nbStates,,drop=FALSE],1,function(x) which(x!=0))))
                 sdind<-unique(unlist(apply(dm[nbStates+1:nbStates,,drop=FALSE],1,function(x) which(x!=0))))
                 if(any(intersect(meanind,sdind))) stop("'DM' for ",i," cannot have parameters in common for mean and concentration")
               }
               tmpbounds <- rbind(matrix(rep(c(-pi,pi),nbStates),ncol=2,byrow=TRUE),matrix(rep(c(0,1),nbStates),ncol=2,byrow=TRUE))
               parNames[[i]] <- c("mean","concentration")
             }
             else {
               parSize[[i]] <- 1
               tmpbounds <- matrix(rep(c(0, 1), parSize[[i]] * nbStates), ncol = 2,byrow = TRUE)
               parNames[[i]] <- c("concentration")
             }
           }
          )
    if(zeroInflation[[i]]) {
      tmpbounds[(parSize[[i]] * nbStates)-nbStates*oneInflation[[i]]-nbStates:1+1,2] <- 1
      parNames[[i]] <- c(parNames[[i]],"zeromass")
    } 
    if(oneInflation[[i]]) {
      tmpbounds[(parSize[[i]] * nbStates)-nbStates:1+1,2] <- 1
      parNames[[i]] <- c(parNames[[i]],"onemass")
    } 
    rownames(tmpbounds) <- paste0(rep(parNames[[i]],each=nbStates),"_",1:nbStates)
    #if(!(dist[[i]] %in% mvndists)) rownames(tmpbounds) <- paste0(rep(parNames[[i]],each=nbStates),"_",1:nbStates)
    #else if(dist[[i]]=="mvnorm2" || dist[[i]]=="rw_mvnorm2"){
    #  rownames(tmpbounds) <- c(paste0(rep(c("mean.x","mean.y"),each=nbStates),"_",1:nbStates),
    #                           paste0(rep(c("sigma.x","sigma.xy","sigma.y"),each=nbStates),"_",1:nbStates))
    #} else if(dist[[i]]=="mvnorm3" || dist[[i]]=="rw_mvnorm3"){
    #  rownames(tmpbounds) <- c(paste0(rep(c("mean.x","mean.y","mean.z"),each=nbStates),"_",1:nbStates),
    #                           paste0(rep(c("sigma.x","sigma.xy","sigma.xz","sigma.y","sigma.yz","sigma.z"),each=nbStates),"_",1:nbStates))
    #}
    bounds[[i]] <- tmpbounds
  }
  
  Bndind <- vector('list',length(distnames))
  names(Bndind) <- distnames
  for(i in distnames){
    Bndind[[i]] <- ifelse(is.null(DM[[i]]),TRUE,FALSE)
  }
  
  if(!is.null(userBounds)) {
    if(!is.list(userBounds) | is.null(names(userBounds))) stop("'userBounds' must be a named list")
    if(!any(names(userBounds) %in% distnames)) stop("userBounds names must include at least one of: ",paste0(distnames,collapse=", "))
    for(i in distnames){
      if(is.null(userBounds[[i]])) 
        userBounds[[i]]<-bounds[[i]]
      else {
        if(!all(dim(bounds[[i]])==dim(userBounds[[i]]))) stop("userBounds for ",i," must be of dimension ",dim(bounds[[i]])[1],"x",dim(bounds[[i]])[2])
        rownames(userBounds[[i]])<-rownames(bounds[[i]])
        Bndind[[i]] <- (isTRUE(all.equal(userBounds[[i]],bounds[[i]])) & Bndind[[i]])
        if(dist[[i]]=="wrpcauchy") bounds[[i]][nrow(bounds[[i]])-(nbStates-1):0,1] <- -1
        if(any(userBounds[[i]][,1]<bounds[[i]][,1]) | any(userBounds[[i]][,2]>bounds[[i]][,2])) stop("userBounds for ",i," must be within the parameter space")
        if(any(userBounds[[i]][,1]>userBounds[[i]][,2])) stop("check userBounds for ",i)
      }
    }
    bounds <- userBounds
  }
  #boundInd <- lapply(DM,getboundInd)
  return(list(parSize=parSize[distnames],bounds=bounds[distnames],parNames=parNames[distnames],Bndind=Bndind[distnames]))
}

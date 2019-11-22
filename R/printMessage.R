printMessage <- function(nbStates,dist,p,DM,formula,formulaDelta,formulaPi,mixtures,what="Fitting",hierarchical=FALSE){
  distnames <- names(dist)
  message("=======================================================================")
  if(!hierarchical) message(what," HMM with ",nbStates," state",ifelse(nbStates>1,"s","")," and ",length(distnames)," data stream",ifelse(length(distnames)>1,"s",""))
  else message(what," hierarchical HMM with ",nbStates," state",ifelse(nbStates>1,"s","")," and ",length(distnames)," data stream",ifelse(length(distnames)>1,"s",""))
  message("-----------------------------------------------------------------------\n")
  for(i in distnames){
    pNames<-p$parNames[[i]]
    #if(!isFALSE(inputs$circularAngleMean[[i]])){ 
    #pNames[1]<-paste0("circular ",pNames[1])
    #if(inputs$consensus[[i]]) pNames[2]<-paste0("consensus ",pNames[2])
    #}
    if(is.null(DM[[i]])){
      message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,"=~1",collapse=", "),")")
    } else if(is.list(DM[[i]])){
      message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,"=",DM[[i]],collapse=", "),")")
    } else message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,": custom",collapse=", "),")")
  }
  message("\n Transition probability matrix formula: ",paste0(formula,collapse=""))
  message("\n Initial distribution formula: ",paste0(formulaDelta,collapse=""))
  if(mixtures>1) {
    message("\n Number of mixtures: ",mixtures)
    message(" Mixture probability formula: ",paste0(formulaPi,collapse=""))
  }
  message("=======================================================================")
}
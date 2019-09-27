
#' Print \code{momentuHMM}
#' @method print momentuHMM
#'
#' @param x A \code{momentuHMM} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' print(m)
#'
#' @export
print.momentuHMM <- function(x,...)
{
  m <- x
  distnames <- names(m$conditions$dist)
  DMind <- m$conditions$DMind
  
  if(!is.null(m$modelName)) {
    mess <- paste("Model:",m$modelName)
    cat(rep("-",nchar(mess)),"------------\n",sep="")
    cat(mess,"\n")
    cat(rep("-",nchar(mess)),"------------\n\n",sep="")
  }
  
  if(!is.null(m$mod$minimum))
    cat("Value of the maximum log-likelihood:",-m$mod$minimum,"\n\n")
  
  for(i in distnames){
    cat("\n")
    if(DMind[[i]]) {
      cat(i,                "parameters:\n")      
      cat(rep("-",nchar(i)),"------------\n",sep="")
      print(m$mle[[i]])
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat(rep("-",nchar(i)),"----------------------------------\n",sep="")
      print(m$CIbeta[[i]]$est)
    }

    if(!DMind[[i]]){
      cat("\n")
      cat(i,                "parameters (based on mean covariate values):\n")
      cat(rep("-",nchar(i)),"---------------------------------------------\n",sep="")
      print(x$CIreal[[i]]$est)
    }
  }
  
  if(length(m$stateNames)>1){
    
    if(!is.null(m$conditions$recharge)){
      cat("\n")
      cat("Recharge parameters for the transition probabilities:\n")
      cat("-----------------------------------------------------\n")
      g0theta <- c(m$mle$g0,m$mle$theta)
      names(g0theta) <- c(paste0("g0:",names(m$mle$g0)),paste0("theta:",names(m$mle$theta)))
      print(g0theta)      
    }
    #if(!is.null(m$mle$beta)) {
      cat("\n")
      cat("Regression coeffs for the transition probabilities:\n")
      cat("---------------------------------------------------\n")
      print(m$mle$beta)
    #}
      
    if(!is.null(m$mle$pi)){
      cat("\n")
      cat("Mixture probabilities:\n")
      cat("----------------------\n")
      if(is.null(m$conditions$formulaPi)) {
        formPi <- ~1
      } else formPi <- m$conditions$formulaPi
      if(!length(attr(terms.formula(formPi),"term.labels")) & is.null(m$conditions$formulaPi)){
        tmp <- m$mle$pi[1,]
        names(tmp) <- paste0("mix",1:m$conditions$mixtures)
        print(tmp)
      } else print(m$mle$pi)
    }
  
    if(!is.null(m$mle$gamma)) {
      cat("\n")
      cat("Transition probability matrix:\n")
      cat("------------------------------\n")
      print(m$mle$gamma)
    } else {
      cat("\n")
      cat("Transition probability matrix (based on mean covariate values):\n")
      cat("---------------------------------------------------------------\n")
      print(m$CIreal$gamma$est)    
    }
  
    cat("\n")
    cat("Initial distribution:\n")
    cat("---------------------\n")
    m <- delta_bc(m)
    if(is.null(m$conditions$formulaDelta)) {
      formDelta <- ~1
    } else formDelta <- m$conditions$formulaDelta
    if(!length(attr(terms.formula(formDelta),"term.labels")) & is.null(m$conditions$formulaDelta)){
      tmp <- m$mle$delta[seq(1,nrow(m$mle$delta),nrow(m$mle$delta)/m$conditions$mixtures),]
      if(m$conditions$mixtures==1) rownames(tmp)<-NULL
      else rownames(tmp) <- paste0("mix",1:m$conditions$mixtures)
      print(tmp)
    } else print(m$mle$delta)
  }
}

#' @method print momentuHierHMM
#' @rdname print.momentuHMM
#' @export
print.momentuHierHMM <- function(x,...)
{
  m <- x
  distnames <- names(m$conditions$dist)
  DMind <- m$conditions$DMind
  
  if(!is.null(m$modelName)) {
    mess <- paste("Model:",m$modelName)
    cat(rep("-",nchar(mess)),"------------\n",sep="")
    cat(mess,"\n")
    cat(rep("-",nchar(mess)),"------------\n\n",sep="")
  }
  
  if(!is.null(m$mod$minimum))
    cat("Value of the maximum log-likelihood:",-m$mod$minimum,"\n\n")
  
  for(i in distnames){
    cat("\n")
    if(DMind[[i]]) {
      cat(i,                "parameters:\n")      
      cat(rep("-",nchar(i)),"------------\n",sep="")
      print(m$mle[[i]])
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat(rep("-",nchar(i)),"----------------------------------\n",sep="")
      print(m$CIbeta[[i]]$est)
    }
    
    if(!DMind[[i]]){
      cat("\n")
      cat(i,                "parameters (based on mean covariate values):\n")
      cat(rep("-",nchar(i)),"---------------------------------------------\n",sep="")
      print(x$CIreal[[i]]$est)
    }
  }
  
  if(length(m$stateNames)>1){
    
    if(!is.null(m$mle$pi)){
      cat("\n")
      cat("Mixture probabilities:\n")
      cat("----------------------\n")
      if(is.null(m$conditions$formulaPi)) {
        formPi <- ~1
      } else formPi <- m$conditions$formulaPi
      if(!length(attr(terms.formula(formPi),"term.labels")) & is.null(m$conditions$formulaPi)){
        tmp <- m$mle$pi[1,]
        names(tmp) <- paste0("mix",1:m$conditions$mixtures)
        print(tmp)
      } else print(m$mle$pi)
    }
    
    hierStates <- m$conditions$hierStates
    hierBeta <- m$conditions$hierBeta
    
    if(!is.list(hierBeta)){
      beta0 <- list(beta=hierBeta)
    } else {
      beta0 <- hierBeta
    }
    delta0 <- m$conditions$hierDelta
    
    if(!is.null(m$conditions$recharge)){
      g0 <- m$mle$g0
      theta <- m$mle$theta
      cat("\n\n")
      cat("---------------------------------------------------------------\n")
      cat("Initial recharge parameter (g0):\n")
      cat("---------------------------------------------------------------\n")
      for(j in 1:(hierStates$height-1)){
        tmpPar <- g0[names(beta0$beta[[paste0("level",j)]]$g0)]
        if(length(tmpPar)){
          cat("-------------------------- ",paste0("level",j)," ---------------------------\n")
          print(tmpPar) 
        }
      }
      cat("---------------------------------------------------------------\n")
      cat("\n")
      cat("---------------------------------------------------------------\n")
      cat("Recharge function parameters (theta):\n")
      cat("---------------------------------------------------------------\n")
      for(j in 1:(hierStates$height-1)){
        tmpPar <- theta[names(beta0$beta[[paste0("level",j)]]$theta)]
        if(length(tmpPar)){
          cat("-------------------------- ",paste0("level",j)," ---------------------------\n")
          print(tmpPar) 
        }
      }
      cat("---------------------------------------------------------------\n")
    } else cat("\n")
    
    cat("\n")
    cat("---------------------------------------------------------------\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("---------------------------------------------------------------\n")
    
    for(j in 1:(hierStates$height-1)){
      cat("-------------------------- ",paste0("level",j)," ---------------------------\n")
      if(j>1){
        for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
          tmpPar <- beta0$beta[[paste0("level",j)]][[jj]]$beta
          print(tmpPar)
          cat("\n")
        }
      } else {
        tmpPar <- beta0$beta[[paste0("level",j)]]$beta
        print(tmpPar)
        cat("\n")
      }
    }
    cat("---------------------------------------------------------------\n")
    
    cat("\n")
    cat("---------------------------------------------------------------\n")
    cat("Transition probability matrix (based on mean covariate values):\n")
    cat("---------------------------------------------------------------\n")
    for(j in 1:(hierStates$height-1)){
      cat("-------------------------- ",paste0("level",j)," ---------------------------\n")
      if(j>1){
        for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
          tmpPar <- m$CIreal$hierGamma[[paste0("level",j)]]$gamma[[jj]]$est
          print(tmpPar)
          cat("\n")
        }
      } else {
        tmpPar <- m$CIreal$hierGamma[[paste0("level",j)]]$gamma$est
        print(tmpPar)
        cat("\n")
      }
    } 
    cat("---------------------------------------------------------------\n")
    
    cat("\n")
    cat("--------------------------------------------------\n")
    cat("Regression coeffs for the initial distribution:\n")
    cat("--------------------------------------------------\n")
    for(j in 1:(hierStates$height-1)){
      cat("-------------------- ",paste0("level",j)," --------------------\n")
      if(j>1){
        for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
          tmpPar <- delta0[[paste0("level",j)]][[jj]]$delta
          print(tmpPar)
          cat("\n")
        }
      } else {
        tmpPar <- delta0[[paste0("level",j)]]$delta
        print(tmpPar)
        cat("\n")
      }
    }
    cat("--------------------------------------------------\n")
    
    cat("\n")
    cat("--------------------------------------------------\n")
    cat("Initial distribution:\n")
    cat("--------------------------------------------------\n")
    for(j in 1:(hierStates$height-1)){
      cat("-------------------- ",paste0("level",j)," --------------------\n")
      if(j>1){
        for(jj in hierStates$Get("name",filterFun=function(x) x$level==j & x$count>0)){
          tmpPar <- m$CIreal$hierDelta[[paste0("level",j)]]$delta[[jj]]$est
          print(tmpPar)
          cat("\n")
        }
      } else {
        tmpPar <- m$CIreal$hierDelta[[paste0("level",j)]]$delta$est
        print(tmpPar)
        cat("\n")
      }
    } 
    cat("--------------------------------------------------\n")
  }
}

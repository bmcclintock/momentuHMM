
context("w2n")

test_that("The output is in the right format",{
  
  m<-example$m
  nbStates <- 2
  nbCovs <- 2
  parSize <- list(step=2,angle=2)
  par <- list(step=c(t(m$mle$step)),angle=c(t(m$mle$angle)))
  bounds <- m$conditions$bounds
  beta <- matrix(rnorm(6),ncol=2,nrow=3)
  delta <- c(0.6,0.4)
  
  distnames<-names(m$conditions$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
  }
  
  wpar <- n2w(par,bounds,beta,delta,nbStates,m$conditions$estAngleMean,NULL,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)
  p <-   w2n(wpar,bounds,parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,1,m$conditions$dist,m$conditions$Bndind,nc,meanind)
  
  expect_equal(length(p$step),parSize$step*nbStates)
  expect_equal(length(p$angle),parSize$angle*nbStates)
  expect_equal(dim(p$beta),dim(beta))
  expect_equal(length(p$delta),length(delta))
})

test_that("w2n and n2w are inverse",{
  m<-example$m
  nbStates <- 2
  nbCovs <- 2
  parSize <- list(step=2,angle=2)
  par <- list(step=c(t(m$mle$step)),angle=c(t(m$mle$angle)))
  bounds <- m$conditions$bounds
  beta <- matrix(rnorm(6),ncol=2,nrow=3)
  delta <- c(0.6,0.4)
  
  distnames<-names(m$conditions$dist)
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(m$conditions$fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(m$conditions$circularAngleMean[[i]]) meanind[[i]] <- which((apply(m$conditions$fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
  }
  
  wpar <- n2w(par,bounds,beta,delta,nbStates,m$conditions$estAngleMean,NULL,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)
  p <-   w2n(wpar,bounds,parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,m$conditions$DMind,m$conditions$workcons,1,m$conditions$dist,m$conditions$Bndind,nc,meanind)
  

  expect_equal(p$step[,1],par$step,tolerance=1e-10)
  expect_equal(p$angle[,1],par$angle,tolerance=1e-10)
  expect_equal(p$beta,beta,tolerance=1e-10)
  expect_equal(p$delta,delta,tolerance=1e-10)
})


context("n2w")

test_that("Exception is thrown",{
  
  m<-example$m
  nbStates <- 2
  par <- list(step=c(0.5,1.5,10,100))
  bounds <- list(step=matrix(c(0,1,0,1,
                     0,Inf,0,Inf),
                   byrow=TRUE,ncol=2))
  beta <- matrix(rnorm(6),ncol=2,nrow=3)
  delta <- c(0.6,0.4)
  
  expect_that(n2w(par,bounds,beta,delta,nbStates,list(step=FALSE),NULL,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind),
              throws_error())
})

test_that("Lengths of input and output are the same",{

  m<-example$m
  nbStates <- 2
  par <- list(step=c(0.5,0.2,10,100))
  bounds <- list(step=matrix(c(0,1,0,1,
                               0,Inf,0,Inf),
                             byrow=TRUE,ncol=2))
  beta <- matrix(rnorm(6),ncol=2,nrow=3)
  delta <- c(0.6,0.4)

  expect_equal(length(n2w(par,bounds,beta,log(delta[-1]/delta[1]),nbStates,list(step=FALSE),NULL,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)),
               length(par$step)+length(beta)+length(delta)-1)
})


context("n2w")

test_that("Exception is thrown",{
  nbStates <- 2
  par <- c(0.5,1.5,10,100)
  bounds <- matrix(c(0,1,0,1,
                     0,Inf,0,Inf),
                   byrow=TRUE,ncol=2)
  beta <- matrix(rnorm(6),ncol=2,nrow=3)
  delta <- c(0.6,0.4)

  expect_that(n2w(par,bounds,beta,delta,nbStates,estAngleMean=FALSE),
              throws_error("Check the parameters bounds."))
})

test_that("Lengths of input and output are the same",{
  nbStates <- 2
  par <- c(0.5,0.2,10,100)
  bounds <- matrix(c(0,1,0,1,
                     0,Inf,0,Inf),
                   byrow=TRUE,ncol=2)
  beta <- matrix(rnorm(6),ncol=2,nrow=3)
  delta <- c(0.6,0.4)

  expect_equal(length(n2w(par,bounds,beta,delta,nbStates,estAngleMean=FALSE)),
               length(par)+length(beta)+length(delta)-1)
})


context("momentuHMM")

test_that("Exceptions are thrown",{
  data <- data.frame(1,2,3)
  states <- rep(1,10)
  mle <- list(1,2,3)
  stepDist <- "gamma"
  angleDist <- "vm"
  mod <- list(1,2,3)
  conditions <- list(dist=list(step=stepDist,angle=angleDist),zeroInflation=list(step=FALSE),estAngleMean=list(angle=TRUE),stationary=FALSE,formula=~c1)
  rawCovs <- data.frame(c1=rnorm(10))

  m <- list(data=data,mle=mle,mod=mod,
            conditions=conditions,rawCovs=rawCovs)
  expect_that(momentuHMM(m),not(throws_error()))

  m <- list(data=data,states=states,mle=mle,stepDist=stepDist)
  expect_that(momentuHMM(m),throws_error())
})

test_that("The output has the right class attribute",{
  data <- data.frame(1,2,3)
  states <- rep(1,10)
  mle <- list(1,2,3)
  stepDist <- "gamma"
  angleDist <- "vm"
  mod <- list(1,2,3)
  conditions <- list(dist=list(step=stepDist,angle=angleDist),zeroInflation=list(step=FALSE),estAngleMean=list(angle=TRUE),stationary=FALSE,formula=~c1)
  rawCovs <- data.frame(c1=rnorm(10))

  m <- list(data=data,mle=mle,mod=mod,
            conditions=conditions,rawCovs=rawCovs)
  m <- momentuHMM(m)

  expect_equal(length(which(class(m)=="momentuHMM")),1)
})

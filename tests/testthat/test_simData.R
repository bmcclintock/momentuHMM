
context("simData")

test_that("Exceptions are thrown",{
  stepPar <- c(1,10,1,5,0.2,0.3)
  anglePar <- c(0,pi,0.5,2)

  expect_that(simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              not(throws_error()))
  expect_that(simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar[1:4],angle=anglePar),nbCovs=2),
              not(throws_error()))

  expect_that(simData(0,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error("nbAnimals should be at least 1."))
  expect_that(simData(1,0,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error("nbStates should be at least 1."))

  expect_that(simData(1,2,dist=list(step="norm",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())
  expect_that(simData(1,2,dist=list(step="norm",angle="norm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())
  expect_that(simData(1,2,dist=list(step="vm",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())
  expect_that(simData(1,2,dist=list(step="gamma",angle="weibull"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())

  stepPar2 <- c(stepPar,1)
  expect_that(simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar2,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())
  anglePar2 <- c(anglePar,1)
  expect_that(simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar2),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())

  stepPar <- c(-1,10,1,5,0.2,0.3)
  expect_that(simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())
  stepPar <- c(1,10,1,-0.5,0.2,0.3)
  expect_that(simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())
  stepPar <- c(1,10,1,5,4,0.3)
  expect_that(simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE)),
              throws_error())
})

test_that("The right slots are defined",{
  stepPar <- c(1,10,1,5,0.2,0.3)
  anglePar <- c(0,pi,0.5,2)
  nbCovs <- 2
  data <- simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE))

  expect_that(!is.null(data$ID),is_true())
  expect_that(!is.null(data$x),is_true())
  expect_that(!is.null(data$y),is_true())
  expect_that(!is.null(data$step),is_true())
  expect_that(!is.null(data$angle),is_true())

  expect_equal(length(which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                              names(data)!="step" & names(data)!="angle")),nbCovs)
})

test_that("The returned object is of the correct class",{
  stepPar <- c(1,10,1,5,0.2,0.3)
  anglePar <- c(0,pi,0.5,2)
  data <- simData(1,2,dist=list(step="gamma",angle="vm"),Par=list(step=stepPar,angle=anglePar),nbCovs=2,zeroInflation=list(step=TRUE))

  expect_equal(class(data),c("momentuHMMData","data.frame"))
})

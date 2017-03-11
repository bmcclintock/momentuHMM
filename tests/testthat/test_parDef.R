
context("parDef")

test_that("Exceptions are thrown",{
  nbStates <- 2
  
  expect_error(parDef(list(step="gamma",angle="vm"),nbStates,list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),NULL,NULL),NA)

  expect_error(parDef(list(step="unif",angle="vm"),nbStates,list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),NULL,NULL))
  expect_error(parDef(list(step="gamma",angle="norm"),nbStates,list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),NULL,NULL))
})

test_that("The output has the right format",{
  nbStates <- 2
  p <- parDef(list(step="gamma",angle="vm"),nbStates,list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),NULL,NULL)
  expect_equal(length(p$parSize),2)
  expect_equal(lapply(p$bounds,nrow),lapply(p$parSize,function(x) x*nbStates))
  expect_equal(lapply(p$bounds,ncol),list(step=2,angle=2))
  expect_equal(lapply(p$parNames,length),p$parSize)

  nbStates <- 3
  p <- parDef(list(step="exp",angle="wrpcauchy"),nbStates,list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),list(step=FALSE,angle=FALSE),NULL,NULL)
  expect_equal(length(p$parSize),2)
  expect_equal(lapply(p$bounds,nrow),lapply(p$parSize,function(x) x*nbStates))
  expect_equal(lapply(p$bounds,ncol),list(step=2,angle=2))
  expect_equal(lapply(p$parNames,length),p$parSize)
})

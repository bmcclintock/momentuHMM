
context("parDef")

test_that("Exceptions are thrown",{
  nbStates <- 2
  expect_that(parDef("gamma","vm",nbStates,FALSE,FALSE),not(throws_error()))

  expect_that(parDef("unif","vm",nbStates,FALSE,FALSE),throws_error())
  expect_that(parDef("gamma","norm",nbStates,FALSE,FALSE),throws_error())
})

test_that("The output has the right format",{
  nbStates <- 2
  p <- parDef("gamma","vm",nbStates,FALSE,FALSE)
  expect_equal(length(p$parSize),2)
  expect_equal(nrow(p$bounds),sum(p$parSize)*nbStates)
  expect_equal(ncol(p$bounds),2)
  expect_equal(length(p$parNames),p$parSize[1])

  nbStates <- 3
  p <- parDef("exp","wrpcauchy",nbStates,FALSE,FALSE)
  expect_equal(length(p$parSize),2)
  expect_equal(nrow(p$bounds),sum(p$parSize)*nbStates)
  expect_equal(ncol(p$bounds),2)
  expect_equal(length(p$parNames),p$parSize[1])
})

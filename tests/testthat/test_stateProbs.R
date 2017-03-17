
context("stateProbs")

test_that("Exceptions are thrown",{
  m <- example$m
  expect_error(stateProbs(m),NA)
})

test_that("State probabilities sum to one",{
  m <- example$m
  sP <- stateProbs(m)
  expect_equal(rowSums(sP),rep(1,nrow(m$data)))
})
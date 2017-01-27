
context("viterbi")

test_that("Output has the right format",{
  m <- example$m
  expect_equal(length(viterbi(m)),nrow(m$data))
})

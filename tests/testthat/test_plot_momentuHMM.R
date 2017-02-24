
context("plot.momentuHMM")

test_that("Exceptions are thrown",{
  expect_error(plot.momentuHMM(example$m,ask=FALSE,plotCI=TRUE),NA)
})
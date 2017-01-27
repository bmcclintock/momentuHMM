
context("logAlpha")

test_that("Exceptions are thrown",{
  m <- example$m
  expect_error(logAlpha(m),NA)

  expect_error(logAlpha(1))
})

test_that("Output has the right format",{
  m <- example$m
  la <- logAlpha(m)

  expect_equal(nrow(la),nrow(m$data))
  expect_equal(ncol(la),ncol(m$mle$step))
})

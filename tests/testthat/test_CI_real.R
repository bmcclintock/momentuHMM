
context("CI_real")

test_that("Output has the right format",{
  m <- example$m
  c <- CI_real(m)

  expect_equal(length(c),4)
  expect_equal(dim(c$step$est),dim(m$mle$step))
  expect_equal(dim(c$step$se),dim(m$mle$step))
  expect_equal(dim(c$step$lower),dim(m$mle$step))
  expect_equal(dim(c$step$upper),dim(m$mle$step))
  expect_equal(dim(c$angle$est),dim(m$mle$angle))
  expect_equal(dim(c$angle$se),dim(m$mle$angle))
  expect_equal(dim(c$angle$lower),dim(m$mle$angle))
  expect_equal(dim(c$angle$upper),dim(m$mle$angle))
})

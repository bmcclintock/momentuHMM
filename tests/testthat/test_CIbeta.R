
context("CIbeta")

test_that("Output has the right format",{
  m <- example$m
  c <- CIbeta(m)

  expect_equal(length(c),4)
  expect_equal(length(c$step$est),length(m$mle$step))
  expect_equal(length(c$step$se),length(m$mle$step))
  expect_equal(length(c$step$lower),length(m$mle$step))
  expect_equal(length(c$step$upper),length(m$mle$step))
  expect_equal(length(c$angle$est),length(m$mle$angle))
  expect_equal(length(c$angle$se),length(m$mle$angle))
  expect_equal(length(c$angle$lower),length(m$mle$angle))
  expect_equal(length(c$angle$upper),length(m$mle$angle))
  expect_equal(length(c$beta$est),length(m$mle$beta))
  expect_equal(length(c$beta$se),length(m$mle$beta))
  expect_equal(length(c$beta$lower),length(m$mle$beta))
  expect_equal(length(c$beta$upper),length(m$mle$beta))
})

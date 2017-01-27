
context("CI")

test_that("Output has the right format",{
  m <- example$m
  c <- CI(m)

  expect_equal(length(c),3)
  expect_equal(dim(c$stepPar$lower),dim(m$mle$stepPar))
  expect_equal(dim(c$stepPar$upper),dim(m$mle$stepPar))
  expect_equal(dim(c$anglePar$lower),dim(m$mle$anglePar))
  expect_equal(dim(c$anglePar$upper),dim(m$mle$anglePar))
  expect_equal(dim(c$beta$lower),dim(m$mle$beta))
  expect_equal(dim(c$beta$upper),dim(m$mle$beta))
})

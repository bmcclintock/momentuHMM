
context("CIreal")

test_that("Output has the right format",{
  m <- example$m
  c <- CIreal(m)

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

test_that("'covs' working correctly",{
  m <- example$m
  
  c1 <- CIreal(m)
  c2 <- CIreal(m,covs=data.frame(cov1=mean(m$data$cov1),cov2=mean(m$data$cov2)))
  expect_equal(c1,c2)
  
  c3 <- CIreal(m,covs=data.frame(cov2=mean(m$data$cov2)))
  expect_equal(c1,c3)
})

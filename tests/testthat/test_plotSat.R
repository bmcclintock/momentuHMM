
#context("plotSat")

#test_that("Exceptions are thrown",{
#  expect_error(plotSat(data.frame(x=rnorm(10),y=rnorm(10)),zoom=3,location=c(0,0),ask=FALSE),NA)
#})
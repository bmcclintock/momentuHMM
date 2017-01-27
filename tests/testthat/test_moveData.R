
context("momentuHMMData")

test_that("Exceptions are thrown",{
  ID <- rep("Animal1",10)
  x <- rep(1,10)
  y <- rep(1,10)
  step <- rep(1,10)

  data <- data.frame(ID=ID,x=x,y=y)
  expect_that(momentuHMMData(data),throws_error())

  data <- cbind(data,step)
  expect_that(momentuHMMData(data),not(throws_error()))
})

test_that("The output has the right class attribute",{
  ID <- rep("Animal1",10)
  x <- rep(1,10)
  y <- rep(1,10)
  step <- rep(1,10)
  data <- data.frame(ID=ID,x=x,y=y,step=step)

  md <- momentuHMMData(data)

  expect_equal(length(which(class(md)=="momentuHMMData")),1)
})

context("moveData")

test_that("Exceptions are thrown",{
  ID <- rep("Animal1",10)
  x <- rep(1,10)
  y <- rep(1,10)
  step <- rep(1,10)
  
  data <- data.frame(ID=ID,x=x,y=y)
  expect_error(moveData(data))
  
  data <- cbind(data,step)
  expect_error(moveData(data),NA)
})

test_that("The output has the right class attribute",{
  ID <- rep("Animal1",10)
  x <- rep(1,10)
  y <- rep(1,10)
  step <- rep(1,10)
  data <- data.frame(ID=ID,x=x,y=y,step=step)
  
  md <- moveData(data)
  
  expect_equal(length(which(class(md)=="moveData")),1)
})
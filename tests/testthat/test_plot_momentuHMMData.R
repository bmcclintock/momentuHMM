
context("plot.momentuHMMData")

test_that("Exceptions are thrown",{
   ID <- rep("Animal1",100)
   x <- seq(-5,5,length=100)+rnorm(100)
   y <- seq(10,12,length=100)+rnorm(100)
   step <- rgamma(100,2,0.25)
   angle <- rvm(100,0,1.5)-pi
   data <- data.frame(ID=ID,x=x,y=y,step=step,angle=angle)

   expect_error(plot.momentuHMMData(data,ask=FALSE),NA)

   data <- data.frame()
   expect_error(plot.momentuHMMData(data,ask=FALSE),"The data input is empty.")

   data <- data.frame(x=x,y=y)
   expect_error(plot.momentuHMMData(data,ask=FALSE))
})

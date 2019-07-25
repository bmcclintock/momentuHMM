
context("examples")

test_that("all examples work correctly",{

  expect_error(devtools::run_examples(run=TRUE),NA)
  
  unlink(paste0("Rplots",1:47,".pdf"))
})
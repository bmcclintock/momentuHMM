context("examples")

test_that("all examples work correctly",{
  
  expect_error(devtools::run_examples(run_dontrun = TRUE,run_donttest = TRUE),NA)
  
  unlink(paste0("Rplots",1:46,".pdf"))
})
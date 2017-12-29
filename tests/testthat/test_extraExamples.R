
context("extraExamples")

test_that("all examples in extraExamples.R work correctly",{
  
  expect_error(source("extraExamples.R"),NA)
  
  load("extraExamples.RData")
  
  mods <- names(oldExample)
  
  for(i in mods){
    for(j in 1:length(oldExample[[i]])){
      expect_equal(abs(oldExample[[i]][[j]]-newExample[[i]][[j]])<1.e-6,TRUE,info=paste(i,ifelse(i=="miFits",j,"")))
    }
  }
  
})
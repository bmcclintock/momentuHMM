
context("extraExamples")

test_that("all examples in extraExamples.R work correctly",{
  
  expect_error(source("extraExamples.R"),NA)
  
  load("extraExamples.RData")
  
  mods <- c(paste0("mod",1:7),"miFits","bestFit","bestFit2")
  
  for(i in mods){
    if(i=="miFits"){
      for(j in 1:length(example[[i]]$HMMfits)){
        expect_equal(abs(example[[i]]$HMMfits[[j]]$mod$minimum-newExample[[i]]$HMMfits[[j]]$mod$minimum)<1.e-6,TRUE,info=paste(i,j))
      }
    } else {
      expect_equal(abs(example[[i]]$mod$minimum-newExample[[i]]$mod$minimum)<1.e-6,TRUE,info=i)
    }
  }
  
})
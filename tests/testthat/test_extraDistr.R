context("extraDistr")

test_that("extraDistr functions are found",{
  set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
  expect_error(sim<-simData(nbStates=2,nbAnimals=2,obsPerAnimal = 30,dist=list(cat="cat2"),Par=list(cat=c(0.1,0.9))),NA)
  expect_error(fit <- fitHMM(sim,nbStates=2,dist=list(cat="cat2"),Par0=list(cat=c(0.1,0.9))),NA)
  expect_error(pr <- pseudoRes(fit),NA)
  expect_error(st <- viterbi(fit),NA)
  expect_error(plot(fit,ask=FALSE),NA)
})
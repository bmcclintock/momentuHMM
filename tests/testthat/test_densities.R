
context("densities")

test_that("C++ and R density functions are identical",{
  x <- seq(0,100,length=1000)
  # dgamma (note the conversion between mean/sd and shape/rate)
  expect_equal(dgamma(x,0.05^2/0.01^2,0.05/0.01^2),c(dgamma_rcpp(x,matrix(0.05,1,1000),matrix(0.01,1,1000))),tolerance=1e-10)
  expect_equal(dgamma(x,1000^2/1000^2,1000/1000^2),c(dgamma_rcpp(x,matrix(1000,1,1000),matrix(1000,1,1000))),tolerance=1e-10)
  expect_equal(dgamma(x,50^2/10^2,50/10^2),c(dgamma_rcpp(x,matrix(50,1,1000),matrix(10,1,1000))),tolerance=1e-10)

  # dweibull
  expect_equal(dweibull(x,0.05,0.05),c(dweibull_rcpp(x,matrix(0.05,1,1000),matrix(0.05,1,1000))),tolerance=1e-10)
  expect_equal(dweibull(x,10,10),c(dweibull_rcpp(x,matrix(10,1,1000),matrix(10,1,1000))),tolerance=1e-10)
  expect_equal(dweibull(x,1,1.5),c(dweibull_rcpp(x,matrix(1,1,1000),matrix(1.5,1,1000))),tolerance=1e-10)

  # dlnorm
  expect_equal(dlnorm(x,-1000,0.05),c(dlnorm_rcpp(x,matrix(-1000,1,1000),matrix(0.05,1,1000))),tolerance=1e-10)
  expect_equal(dlnorm(x,100,100),c(dlnorm_rcpp(x,matrix(100,1,1000),matrix(100,1,1000))),tolerance=1e-10)
  expect_equal(dlnorm(x,0,2),c(dlnorm_rcpp(x,matrix(0,1,1000),matrix(2,1,1000))),tolerance=1e-10)

  # dexp
  expect_equal(dexp(x,0.05),c(dexp_rcpp(x,matrix(0.05,1,1000),matrix())),tolerance=1e-10)
  expect_equal(dexp(x,100),c(dexp_rcpp(x,matrix(100,1,1000),matrix())),tolerance=1e-10)
  expect_equal(dexp(x,1),c(dexp_rcpp(x,matrix(1,1,1000),matrix())),tolerance=1e-10)
  
  # dpois
  x <- seq(0,999)
  expect_equal(dpois(x,0.05),c(dpois_rcpp(x,matrix(0.05,1,1000),matrix())),tolerance=1e-10)
  expect_equal(dpois(x,100),c(dpois_rcpp(x,matrix(100,1,1000),matrix())),tolerance=1e-10)
  expect_equal(dpois(x,1),c(dpois_rcpp(x,matrix(1,1,1000),matrix())),tolerance=1e-10)
  
  #dbeta
  x <- seq(0,1,length=1000)
  expect_equal(dbeta(x,0.05,0.05),c(dbeta_rcpp(x,matrix(0.05,1,1000),matrix(0.05,1,1000))),tolerance=1e-10)
  expect_equal(dbeta(x,100,100),c(dbeta_rcpp(x,matrix(100,1,1000),matrix(100,1,1000))),tolerance=1e-10)
  expect_equal(dbeta(x,1,1),c(dbeta_rcpp(x,matrix(1,1,1000),matrix(1,1,1000))),tolerance=1e-10)
  
  x <- seq(-pi,pi,length=1000)
  # dvm
  expect_equal(dvm(x,0,0.05),c(dvm_rcpp(x,matrix(0,1,1000),matrix(0.05,1,1000))),tolerance=1e-10)
  expect_equal(dvm(x,pi,1000),c(dvm_rcpp(x,matrix(pi,1,1000),matrix(1000,1,1000))),tolerance=1e-10)
  expect_equal(dvm(x,pi/2,1),c(dvm_rcpp(x,matrix(pi/2,1,1000),matrix(1,1,1000))),tolerance=1e-10)

  # dwrpcauchy
  expect_equal(dwrpcauchy(x,0,0.01),c(dwrpcauchy_rcpp(x,matrix(0,1,1000),matrix(0.01,1,1000))),tolerance=1e-10)
  expect_equal(dwrpcauchy(x,pi,0.99),c(dwrpcauchy_rcpp(x,matrix(pi,1,1000),matrix(0.99,1,1000))),tolerance=1e-10)
  expect_equal(dwrpcauchy(x,pi/2,0.5),c(dwrpcauchy_rcpp(x,matrix(pi/2,1,1000),matrix(0.5,1,1000))),tolerance=1e-10)
  
  
  # dnorm
  x <- seq(-100,100,length=1000)
  expect_equal(dnorm(x,0,1),c(dnorm_rcpp(x,matrix(0,1,1000),matrix(1,1,1000))),tolerance=1e-10)
  expect_equal(dnorm(x,-100,100),c(dnorm_rcpp(x,matrix(-100,1,1000),matrix(100,1,1000))),tolerance=1e-10)
  expect_equal(dnorm(x,500,1000),c(dnorm_rcpp(x,matrix(500,1,1000),matrix(1000,1,1000))),tolerance=1e-10)
  
  #dbern
  x <- c(0,1)
  expect_equal(dbern(x,0.5),c(dbern_rcpp(x,matrix(0.5,1,2),matrix(0,1,2))),tolerance=1e-10)
  expect_equal(dbern(x,0.01),c(dbern_rcpp(x,matrix(0.01,1,2),matrix(0,1,2))),tolerance=1e-10)
  expect_equal(dbern(x,0.95),c(dbern_rcpp(x,matrix(0.95,1,2),matrix(0,1,2))),tolerance=1e-10)
  
  #mvnorm2
  x <- matrix(c(seq(0,1000,length=1000),seq(-100,0,length=1000)),1000,2)
  expect_equal(mvtnorm::dmvnorm(x,sigma=matrix(c(500,0,0,500),2,2)),c(dmvnorm2(x,matrix(0,2,1000),matrix(sqrt(c(500,500,0)),3,2*1000))),tolerance=1e-10)
  expect_equal(mvtnorm::dmvnorm(x,mean=c(100,1000),sigma=matrix(c(10,0,0,10),2,2)),c(dmvnorm2(x,matrix(c(100,1000),2,1000),matrix(sqrt(c(10,10,0)),3,2*1000))),tolerance=1e-10)
  expect_equal(mvtnorm::dmvnorm(x,mean=c(-100,100),sigma=matrix(c(100,-5,-5,100),2,2)),c(dmvnorm2(x,matrix(c(-100,100),2,1000),matrix(c(10,10,-0.5),3,2*1000))),tolerance=1e-10)
  
  #mvnorm3
  x <- matrix(c(seq(0,1000,length=1000),seq(-100,0,length=1000),seq(0,1000,length=1000)),1000,3)
  expect_equal(mvtnorm::dmvnorm(x,sigma=matrix(c(500,0,0,0,500,0,0,0,500),3,3)),c(dmvnorm3(x,matrix(0,3,1000),matrix(sqrt(c(500,500,500,0,0,0)),6,3*1000))),tolerance=1e-10)
  expect_equal(mvtnorm::dmvnorm(x,mean=c(100,1000,100),sigma=matrix(c(10,0,0,0,10,0,0,0,10),3,3)),c(dmvnorm3(x,matrix(c(100,1000,100),3,1000),matrix(sqrt(c(10,10,10,0,0,0)),6,3*1000))),tolerance=1e-10)
  expect_equal(mvtnorm::dmvnorm(x,mean=c(-100,100,-100),sigma=matrix(c(100,-5,-5,-5,100,-5,-5,-5,100),3,3)),c(dmvnorm3(x,matrix(c(-100,100,-100),3,1000),matrix(c(10,10,10,-0.05,-0.05,-0.05),6,3*1000))),tolerance=1e-10)
  
  #mvnorm3
  x <- matrix(c(seq(0,1000,length=1000),seq(-100,0,length=1000),seq(0,1000,length=1000)),1000,3)
  expect_equal(mvtnorm::dmvnorm(x,sigma=matrix(c(500,0,0,0,500,0,0,0,500),3,3)),c(dmvnorm3(x,matrix(0,3,1000),matrix(sqrt(c(500,500,500,0,0,0)),6,3*1000))),tolerance=1e-10)
  expect_equal(mvtnorm::dmvnorm(x,mean=c(100,1000,100),sigma=matrix(c(10,0,0,0,10,0,0,0,10),3,3)),c(dmvnorm3(x,matrix(c(100,1000,100),3,1000),matrix(sqrt(c(10,10,10,0,0,0)),6,3*1000))),tolerance=1e-10)
  expect_equal(mvtnorm::dmvnorm(x,mean=c(-100,100,-100),sigma=matrix(c(100,-5,-5,-5,100,-5,-5,-5,100),3,3)),c(dmvnorm3(x,matrix(c(-100,100,-100),3,1000),matrix(c(10,10,10,-0.05,-0.05,-0.05),6,3*1000))),tolerance=1e-10)
  
  #dcat
  x <- extraDistr::rcat(1000,prob=rep(1/10,10))
  expect_equal(extraDistr::dcat(x,rep(1/10,10)),c(dcat_rcpp(x,matrix(1/10,10,1000),matrix(0,1,2))),tolerance=1e-10)
  set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
  pr <- matrix(runif(1000*10),1000,10)
  pr <- pr/rowSums(pr)
  expect_equal(extraDistr::dcat(x,pr),c(dcat_rcpp(x,t(pr),matrix(0,1,2))),tolerance=1e-10)
  
  # dnegbinom
  x <- seq(0,999)
  expect_equal(dnbinom(x,size=500,mu=500),c(dnbinom_rcpp(x,matrix(500,1,1000),matrix(500,1,1000))),tolerance=1e-10)
  expect_equal(dnbinom(x,size=1,mu=500),c(dnbinom_rcpp(x,matrix(500,1,1000),matrix(1,1,1000))),tolerance=1e-10)
  expect_equal(dnbinom(x,size=1000,mu=1),c(dnbinom_rcpp(x,matrix(1,1,1000),matrix(1000,1,1000))),tolerance=1e-10)
  
  # dlogis
  x <- seq(-100,100,length=1000)
  expect_equal(dlogis(x,0,1),c(dlogis_rcpp(x,matrix(0,1,1000),matrix(1,1,1000))),tolerance=1e-10)
  expect_equal(dlogis(x,-100,100),c(dlogis_rcpp(x,matrix(-100,1,1000),matrix(100,1,1000))),tolerance=1e-10)
  expect_equal(dlogis(x,500,1000),c(dlogis_rcpp(x,matrix(500,1,1000),matrix(1000,1,1000))),tolerance=1e-10)
  
  # dt
  x <- seq(-100,100,length=1000)
  expect_equal(dt(x,ncp=0,df=1),c(dt(x,matrix(1,1,1000),matrix(0,1,1000))),tolerance=1e-10)
  expect_equal(dt(x,ncp=-100,df=100),c(dt_rcpp(x,matrix(100,1,1000),matrix(-100,1,1000))),tolerance=1e-10)
  expect_equal(dt(x,ncp=500,df=1000),c(dt_rcpp(x,matrix(1000,1,1000),matrix(500,1,1000))),tolerance=1e-10)
  
})

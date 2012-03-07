
require(testthat)
require(biwt)

context("ksolve") 

test_that("check that ksolve is pretty close to original code.", {

  m <- matrix(rnorm(20), ncol=2)
  med.init=covMcd(m)
  r <- 0.2
  p<-2
  n <- dim(m)[2]
  c1<-rejpt.bw(p=2,r)[1]
  b0<-erho.bw(p=2,c1)[1]  
  d <- sqrt(mahalanobis(m,med.init$center,med.init$cov))
  knew <- fastksolve(d,p,c1,b0)
  korig <- ksolve(d,p,c1,b0)

  expect_that( abs(knew-korig) < 0.001 , is_true())


})



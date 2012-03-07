
require(testthat)
require(biwt)

context("fastukey parallel correlation") 


test_that("parallel tukey cor returns same matrix regardless of cores used, small matrix.", {
  
  x <- matrix(rnorm(mean=5, sd=2, n=50), nrow=5)
  res1 <- partukeycor(x)
  res2 <- partukeycor(x, cores=2)
  res3 <- partukeycor(x, cores=3)
  
  expect_that( all(dim(res1) == dim(res2))  , is_true())
  expect_that( all(dim(res2) == dim(res3))  , is_true())

  res1 <- res1[upper.tri(res1)]
  res2 <- res2[upper.tri(res2)]
  expect_that(all(res1 == res2), is_true())
  
})



test_that("parallel tukey cor returns result close to biwt.cor", {
  
  x <- matrix(rnorm(mean=5, sd=2, n=1000), nrow=50)
  res0 <- biwt.cor(x)
  res1 <- partukeycor(x, cores=2)
  res0 <- res0[upper.tri(res0)]
  res1 <- res1[upper.tri(res1)]
  
  expect_that( all(abs(res1-res0) < 0.001) , is_true())


})


test_that("parallel tukey cor returns same matrix regardless of cores used, big matrix.", {
  
  x <- matrix(rnorm(mean=5, sd=2, n=1000), nrow=50)
  res1 <- partukeycor(x)
  res2 <- partukeycor(x, cores=2)
  res3 <- partukeycor(x, cores=3)
  
  expect_that( all(dim(res1) == dim(res2))  , is_true())
  expect_that( all(dim(res2) == dim(res3))  , is_true())

})


library(randomForestSRC)
library(purrr)
library(testthat)

test_that("`emd()` works for matrices", {
  set.seed(1)
  X = matrix(rnorm(3*2, mean=-1),ncol=2) # m obs. for X
  Y = matrix(rnorm(5*2, mean=+1),ncol=2) # n obs. for Y
  expect_snapshot(emd(X, Y))
})

test_that("`emd()` works for matrices", {
  set.seed(1)
  X = matrix(rnorm(3*2, mean=-1),ncol=2) # m obs. for X
  expect_snapshot(emd(X, X))
})

test_that("`emd()` works for data.frames", {
  expect_snapshot(emd(iris[1:75,], iris[76:150,]))
})

test_that("`emd()` works for data.frames 2", {
  expect_snapshot(emd(iris[1:75,], iris[1:75,]))
})


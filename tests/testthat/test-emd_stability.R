library(randomForestSRC)
library(purrr)
library(testthat)

test_that("`emd_stability()` works for matrices", {
  set.seed(1)
  X = matrix(rnorm(3*2, mean=-1),ncol=2) # m obs. for X
  Y = matrix(rnorm(5*2, mean=+1),ncol=2) # n obs. for Y
  expect_snapshot(emd_stability(X, Y))
})

test_that("`emd_stability()` works for matrices 2", {
  set.seed(1)
  X = matrix(rnorm(3*2, mean=-1),ncol=2) # m obs. for X
  expect_snapshot(emd_stability(X, X))
})

test_that("`emd_stability()` works for data.frames", {
  expect_snapshot(emd_stability(iris[1:10,], iris[11:20,]))
})

test_that("`emd_stability()` works for data.frames 2", {
  expect_snapshot(emd_stability(iris[1:10,], iris[1:10,]))
})

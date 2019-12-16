# Change this to source your lasso functions
source("LassoFunctions.R")

# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Do 2 tests for soft-thresholding function below
#################################################
outl1 <- soft(-2, 2)
outc1 <- soft_c(-2, 2)

test1 <- tryCatch(test_that("soft is correct", {
  expect_equal(outl1, outc1)
}), error = function(e) 0)


outl2 <- soft(-6, 4)
outc2 <- soft_c(-6, 4)

test2 <- tryCatch(test_that("soft is correct", {
  expect_equal(outl2, outc2)
}), error = function(e) 0)

# Do 2 tests for lasso objective function below
#################################################

set.seed(1)
p <- 50
n <- 30
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda <- 0.3
beta <- rnorm(p)

outl3 <- lasso(out$Xtilde, out$Ytilde, beta, lambda)
outc3 <- lasso_c(out$Xtilde, out$Ytilde, beta, lambda)

test3 <- tryCatch(test_that("objective is correct", {
  expect_equal(as.numeric(outl3), as.numeric(outc3))
}), error = function(e) 0)


set.seed(100)
p <- 100
n <- 200
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda <- 1.3
beta <- rnorm(p)

outl4 <- lasso(out$Xtilde, out$Ytilde, beta, lambda)
outc4 <- lasso_c(out$Xtilde, out$Ytilde, beta, lambda)

test4 <- tryCatch(test_that("objective is correct", {
  expect_equal(as.numeric(outl4), as.numeric(outc4))
}), error = function(e) 0)

# Do 2 tests for fitLASSOstandardized function below
#################################################

set.seed(60)
n <- 50
p <- 30
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda_max <- max(abs(crossprod(out$Xtilde, out$Ytilde)) / nrow(X))
beta_start <- rep(0.1, ncol(X))

outl5 <- fitLASSOstandardized(out$Xtilde, as.vector(out$Ytilde), lambda = lambda_max, beta_start, eps = 1e-6)
outc5 <- fitLASSOstandardized_c(out$Xtilde, as.vector(out$Ytilde), lambda = lambda_max, beta_start, eps = 1e-6)

test5 <- tryCatch(test_that("lambda max gives 0", {
  expect_equal(as.vector(outl5$beta), as.vector(outc5))
}), error = function(e) 0)


set.seed(100)
n <- 50
p <- 30
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda_max <- max(abs(crossprod(out$Xtilde, out$Ytilde)) / nrow(X))
beta_start <- rep(0.1, ncol(X))

outl6 <- fitLASSOstandardized(out$Xtilde, as.vector(out$Ytilde), lambda = 0, beta_start = beta_start, eps = 1e-6)
outc6 <- fitLASSOstandardized_c(out$Xtilde, as.vector(out$Ytilde), lambda = 0, beta_start = beta_start, eps = 1e-6)

test6 <- tryCatch(test_that("lambda 0 gives least square solution", {
  expect_equal(as.vector(outl6$beta), as.vector(outc6))
}), error = function(e) 0)

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

library(microbenchmark)
set.seed(100)
n <- 50
p <- 30
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda_max <- max(abs(crossprod(out$Xtilde, out$Ytilde)) / nrow(X))
beta_start <- rep(0.1, ncol(X))

microbenchmark(
  fitLASSOstandardized(out$Xtilde, as.vector(out$Ytilde), lambda = 0, beta_start = beta_start, eps = 1e-6),
  fitLASSOstandardized_c(out$Xtilde, as.vector(out$Ytilde), lambda = 0, beta_start = beta_start, eps = 1e-6),
  times = 10
)

# Do 2 tests for fitLASSOstandardized_seq function below
#################################################

set.seed(100)
n <- 50
p <- 30
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda_max <- max(abs(crossprod(out$Xtilde, out$Ytilde)) / nrow(X))
lambda_seq <- seq(lambda_max/2, 0, length.out = 10)


outl7 <- fitLASSOstandardized_seq(out$Xtilde, as.vector(out$Ytilde), lambda_seq = lambda_seq, eps = 1e-6)
outc7 <- fitLASSOstandardized_seq_c(out$Xtilde, as.vector(out$Ytilde), lambda_seq = lambda_seq, eps = 1e-6)

test7 <- tryCatch(test_that("beta for sequence of lambdas", {
  expect_equal(outl7$beta_mat, outc7)
}), error = function(e) 0)

set.seed(200)
n <- 50
p <- 30
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda_max <- max(abs(crossprod(out$Xtilde, out$Ytilde)) / nrow(X))
lambda_seq <- seq(lambda_max, 0, length.out = 10)


outl8 <- fitLASSOstandardized_seq(out$Xtilde, as.vector(out$Ytilde), lambda_seq = lambda_seq, eps = 1e-6)
outc8 <- fitLASSOstandardized_seq_c(out$Xtilde, as.vector(out$Ytilde), lambda_seq = lambda_seq, eps = 1e-6)

test8 <- tryCatch(test_that("beta for sequence of lambdas", {
  expect_equal(outl8$beta_mat, outc8)
}), error = function(e) 0)


# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

set.seed(200)
n <- 50
p <- 30
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
out <- standardizeXY(X, Y)
lambda_max <- max(abs(crossprod(out$Xtilde, out$Ytilde)) / nrow(X))
lambda <- seq(lambda_max, 0, length.out = 10)

microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, as.vector(out$Ytilde), lambda = lambda, eps = 1e-6),
  fitLASSOstandardized_seq_c(out$Xtilde, as.vector(out$Ytilde), lambda = lambda, eps = 1e-6),
  times = 10
)

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler for now
outL <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outL$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outL$lambda_seq),
  times = 10
)

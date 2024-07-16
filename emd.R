library(foreach)
library(T4transport)
library(Matrix)
library(randomForestSRC)
library(dplyr)
library(ggplot2)
library(S7)

m = 3
n = 5
X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y

row_outer_product = function(x, y, fun) {
  if (!getDoParRegistered()) {
    registerDoSEQ()
  }
  ret = foreach(i = seq_len(nrow(x)), .combine = rbind) %dopar% {
    m = matrix(0, nrow = 1, ncol = nrow(y))
    for (j in seq_len(nrow(y))) {
      m[1,j] = fun(x[i,], y[j,])
    }
    m
  }
  ret
}

minkowski = function(x, y, p = 2) {
  row_outer_product(x, y, \(x, y) sum((x - y)^2)^(1/p))
}

cost = minkowski(X, Y)

setGeneric(
  "emd", 
  function(x, y, p, keep_pairs, eps) standardGeneric("emd")
)

setMethod(
  "emd",
  signature(
    x = "Matrix",
    y = "Matrix",
    p = "numeric",
    keep_pairs = "ANY",
    eps = "ANY"
  ),
  function(x, y, p, keep_pairs, eps) {
    cost = minkowski(x, y, p)
    t4 = wassersteinD(cost, p = p)
    if (missing(eps)) {
      eps = 1e-8
    }
    if (!is.numeric(eps)) {
      stop("The eps argument must be the numerical precision.")
    }
    if (missing(keep_pairs)) {
      keep_pairs = FALSE
    }
    if (!is.logical(keep_pairs) || length(keep_pairs) != 1) {
      stop("`keep_pairs` should be a single logical value.")
    }
    if (keep_pairs) {
      matches = Reduce(
        cbind, 
        mat2triplet(t4$plan) |> 
          (\(x) {
            keep = which(abs(x$x) > eps)
            list(
              i = x$i[keep],
              j = x$j[keep],
              x = x$x[keep]
            )
          })()
      )
      # xp = x[matches[, 1], ]
      # yp = y[matches[, 2], ]
      colnames(matches) = NULL
      list(
        # Calculated as sum(matches[,3] * apply((xp - yp)^p, 1, sum))^(1/p).
        dist = t4$distance,
        pairs = matches
      )
    } else {
      list(
        dist = t4$dist,
        pairs = NULL
      )
    }
  }
)

setMethod(
  "emd",
  signature(
    x = "Matrix",
    y = "Matrix",
    p = "missing",
    keep_pairs = "ANY",
    eps = "ANY"
  ),
  function(x, y, p, eps) {
    emd(x, y, 2, eps)
  }
)

setMethod(
  "emd",
  signature(
    x = "matrix",
    y = "matrix",
    p = "numeric",
    keep_pairs = "ANY",
    eps = "ANY"
  ),
  function(x, y, p, eps) {
    emd(as(x, "Matrix"), as(y, "Matrix"), p = p, eps = eps)
  }
)

setMethod(
  "emd",
  signature(
    x = "matrix",
    y = "matrix",
    p = "missing",
    keep_pairs = "ANY",
    eps = "ANY"
  ),
  function(x, y, p, eps) {
    emd(as(x, "Matrix"), as(y, "Matrix"), p = 2, eps = eps)
  }
)

# Helper function to turn a data.frame into a model matrix.
df_to_matrix = function(df) {
  mm = as(model.matrix( ~. - 1, df), "Matrix")
  if (nrow(mm) < nrow(df)) {
    warning("Dropping ", nrow(df) - nrow(mm), " rows.")
  }
  mm
}

setMethod(
  "emd",
  signature(
    x = "data.frame",
    y = "data.frame",
    p = "numeric",
    keep_pairs = "ANY",
    eps = "ANY"
  ),
  function(x, y, p) {
    xm = df_to_matrix(x)
    ym = df_to_matrix(y)
    if (ncol(xm) != ncol(ym)) {
      stop("Number of columns differ after making the model matrix.")
    }
    if (!all(names(xm) == names(ym))) {
      stop("Data frame names don't line up.")
    }
    emd(xm, ym, p)
  }
)

setMethod(
  "emd",
  signature(
    x = "data.frame",
    y = "data.frame",
    p = "missing",
    keep_pairs = "ANY",
    eps = "ANY"
  ),
  function(x, y, p, keep_pairs, eps) {
    emd(x, y, 2)
  }
)

emd(X, Y, 2)
emd(X, Y)

iris1 = iris |> filter(Sepal.Length < 5.8)
iris2 = iris |> filter(Sepal.Length > 5.8)

fit1 = rfsrc(
  Species ~ .,
  data = iris1
)

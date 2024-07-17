
# Helper function apply a function to all row combinations.
#' @importFrom purrr map reduce
row_outer_product = function(x, y, fun) {
  map(
      seq_len(nrow(x)),
      ~ {
        m = matrix(0, nrow = 1, ncol = nrow(y))
        for (j in seq_len(nrow(y))) {
          m[1, j] = fun(x[.x, ], y[j, ])
        }
        m
    }
  ) |> reduce(rbind)
}

# The Minkowski Distance
minkowski = function(x, y, p = 2) {
  row_outer_product(x, y, \(x, y) sum((x - y)^2)^(1/p))
}

setClassUnion("numeric_or_missing", c("missing", "numeric"))
setClassUnion("logical_or_missing", c("missing", "logical"))
setOldClass("rfsrc")

#' @title The Earthmover Distance
#' @description Calculate the earthmover distance between data sets `x` and 
#' `y`.
#' @aliases 
#'  emd,data.frame,data.frame,numeric_or_missing-method
#'  emd,matrix,matrix,numeric_or_missing-method
#'  emd,Matrix,Matrix,numeric_or_missing-method
#' @param x a `matrix`, `Matrix`, or `data.frame`.
#' @param y a `matrix`, `Matrix`, or `data.frame`.
#' @param p an exponent for the order of the distance (default: 2)
#' @return a list with the following elements:
#' * dist - the earthmover distance between `x` and `y`.
#' * pairs - the transportation plan from rows of `x` to rows of `y`
#'   along with the amount of mass to move.
#' @examples
#' X = matrix(rnorm(3*2, mean=-1),ncol=2) # m obs. for X
#' Y = matrix(rnorm(5*2, mean=+1),ncol=2) # n obs. for Y
#' emd(X, Y)
#' @importFrom transport transport
#' @importFrom Matrix Matrix
#' @docType methods
#' @rdname emd-methods
#' @export
setGeneric(
  "emd", 
  function(x, y, p) standardGeneric("emd")
)

# The implementation of emd().
emd_impl = function(x, y, p) {
  if (missing(p)) {
    p = 2
  }
  cost = minkowski(x, y, p)
  matches = transport(
    rep(1/nrow(x), nrow(x)),
    rep(1/nrow(y), nrow(y)),
    cost
  )
  xp = x[matches[, 1], ]
  yp = y[matches[, 2], ]
  colnames(matches) = NULL
  list(
    # Log this for numerical stability?
    dist = sum(matches[,3] * apply((xp - yp)^p, 1, sum))^(1/p),
    pairs = matches
  )
}

setMethod(
  "emd",
  signature(
    x = "Matrix",
    y = "Matrix",
    p = "numeric_or_missing"
  ),
  emd_impl
)

setMethod(
  "emd",
  signature(
    x = "matrix",
    y = "matrix",
    p = "numeric_or_missing"
  ),
  function(x, y, p) {
    if (missing(p)) {
      p = 2
    }
    emd(as(x, "Matrix"), as(y, "Matrix"), p = p)
  }
)

# Helper function to turn a data.frame into a model matrix.
#' @importFrom stats model.matrix
df_to_matrix = function(df) {
  mm = as(model.matrix( ~. - 1, df), "Matrix")
  if (nrow(mm) < nrow(df)) {
    warning("Dropping ", nrow(df) - nrow(mm), " rows.")
  }
  mm
}

create_matrices_from_data_frames = function(x, y) {
  xm = df_to_matrix(x)
  ym = df_to_matrix(y)
  if (ncol(xm) != ncol(ym)) {
    stop("Number of columns differ after making the model matrix.")
  }
  if (!all(names(xm) == names(ym))) {
    stop("Data frame names don't line up.")
  }
  list(xm = xm, ym = ym)  
}

setMethod(
  "emd",
  signature(
    x = "data.frame",
    y = "data.frame",
    p = "numeric_or_missing"
  ),
  function(x, y, p) {
    if (missing(p)) {
      p = 2
    }
    l = create_matrices_from_data_frames(x, y)
    emd(l$xm, l$ym, p)
  }
)

#' @title The Earthmover Distance in a Model Embedding
#' @description Calculate the earthmover distance between data sets `x` and
#' `y` in the embedding defined by `model`.
#' @aliases embedded_emd,data.frame,data.frame,rfsrc,numeric_or_missing-method
#' @param x a `matrix`, `Matrix`, or `data.frame`.
#' @param y a `matrix`, `Matrix`, or `data.frame`.
#' @param model the model that will induce the embedded 
#' representation of `x` and `y`.
#' @param p an exponent for the order of the distance (default: 2)
#' @return a list with the following elements:
#' * dist - the earthmover distance between `x` and `y`.
#' * pairs - the transportation plan from rows of `x` to rows of `y`
#'   along with the amount of mass to move.
#' @examples
#' library(dplyr)
#' library(randomForestSRC)
#'
#' # Subset iris for our x and y.
#' iris1 = iris |> filter(Sepal.Length < 5.8)
#' iris2 = iris |> filter(Sepal.Length >= 5.8)
#'
#' # Create a model.
#' fit1 = rfsrc(
#'   Species ~ .,
#'   data = iris1
#' )
#'
#' # Calculate the embedding distance between the data sets.
#' embedded_emd(iris1, iris2, fit1)
#' @importFrom stats cmdscale
#' @docType methods
#' @rdname embedded_emd-methods
#' @export
setGeneric(
  "embedded_emd",
  function(x, y, model, p) standardGeneric("embedded_emd")
)

#' @importFrom stats predict
embed_samples = function(x, y, model) {
  xy = rbind(x, y)
  dists = predict(model, xy, distance = TRUE)$distance
  cs = cmdscale(d = dists, k = ncol(x)-1)
  xm = cs[seq_len(nrow(x)),]
  ym = cs[(nrow(x) + 1):nrow(xy),]
  list(xm = xm, ym = ym)
}

setMethod(
  "embedded_emd",
  signature(
    x = "data.frame",
    y = "data.frame",
    model = "rfsrc",
    p = "numeric_or_missing"
  ),
  function(x, y, model, p) {
    if (missing(p)) {
      p = 2
    }
    l = embed_samples(x, y, model)
    emd(l$xm, l$ym)
  }
)

#' @title Stability of the Earthmove Distance in Each Sample
#' @description The earthmover distance is the sum of the Minkowski 
#' distances of samples, weighted their transport coefficients. The 
#' `emd_stability()` function evaluates how stable the distance is in each
#' of the samples by iteratively removing each sample, calating the new
#' earthmover distance and subtracting that from the earthmover distance
#' will all samples. The difference between those distances can then be used
#' to evaluate both how dis-simmilar the sample is to other samples as well
#' as how much the distance is dependent on the sample.
#' @aliases
#'   emd_stability,data.frame,data.frame,numeric_or_missing,logical_or_missing-method
#'   emd_stability,matrix,matrix,numeric_or_missing,logical_or_missing-method
#'   emd_stability,Matrix,Matrix,numeric_or_missing,logical_or_missing-method
#' @param x a `matrix`, `Matrix`, or `data.frame`.
#' @param y a `matrix`, `Matrix`, or `data.frame`.
#' @param p an exponent for the order of the distance (default: 2).
#' @param progress Should the progress be shown as the calculation is being
#' performed (default: `interactive()`)? 
#' @return a tibble with the following variable.
#' * var: The variable (either "x" or "y").
#' * row: The row for the corresponding variable.
#' * diff: The difference between the earthmover distance and the earthmover
#'   distance with the row removed.
#' * p_overall: the empirical p-value of the difference compared to all other
#'   differences.
#' * p_within: the empirical p-value of the differences compared to other
#'   difference in the same data set.
#' * p_across: the empirical p-value of the difference compared to differences
#'   in the other data set.
#' @examples
#' X = matrix(rnorm(3*2, mean=-1),ncol=2) # m obs. for X
#' Y = matrix(rnorm(5*2, mean=+1),ncol=2) # n obs. for Y
#' emd_stability(X, Y)
#' @docType methods
#' @rdname emd_stability-methods
#' @export
setGeneric(
  "emd_stability",
  function(x, y, p, progress) standardGeneric("emd_stability")
)

#' @importFrom stats ecdf
two_sided_percentile = function(d) {
  vapply(
    seq_along(d), 
    \(x) {
      r = ecdf(d[-x])(d[x])
      min(r, 1-r)
    },
    NA_real_
  )
}

calc_two_sided_percentile_across = function(d) {
  x_diff = d$diff[d$var == "x"]
  y_diff = d$diff[d$var == "y"]
  xp = ecdf(y_diff)(x_diff)
  xp = vapply(xp, \(x) min(x, 1 - x), NA_real_)
  yp = ecdf(x_diff)(y_diff)
  yp = vapply(yp, \(x) min(x, 1 - x), NA_real_)
  c(xp, yp)
}

#' @importFrom furrr future_map_dbl furrr_options
#' @importFrom tibble tibble
setMethod(
  "emd_stability",
  signature(
    x = "Matrix",
    y = "Matrix",
    p = "numeric_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, progress) {
    if (missing(progress)) {
      progress = interactive()
    }
    ref_dist = emd(x, y, p)$dist
    xs = future_map_dbl(
      seq_len(nrow(x)),
      ~ ref_dist - emd_impl(x[-.x,,drop = FALSE], y, p)$dist,
      .options = furrr_options(seed=TRUE),
      .progress = progress
    )
    ys = future_map_dbl(
      seq_len(nrow(y)),
      ~ ref_dist - emd_impl(x, y[-.x,,drop = FALSE], p)$dist,
      .options=furrr_options(seed=TRUE),
      progress = progress
    )
    ret = tibble(
      var = c(rep("x", length(xs)), rep("y", length(ys))),
      row = c(seq_len(length(xs)), seq_len(length(ys))),
      diff = c(xs, ys)
    )
    ret$p_overall = two_sided_percentile(ret$diff)
    ret$p_within = c(
      two_sided_percentile(ret$diff[ret$var == "x"]),
      two_sided_percentile(ret$diff[ret$var == "y"])
    )
    ret$p_across = calc_two_sided_percentile_across(ret)
    ret
  }
)

#' @importFrom methods as
setMethod(
  "emd_stability",
  signature(
    x = "matrix",
    y = "matrix",
    p = "numeric_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(progress)) {
      progress = interactive()
    }
    emd_stability(as(x, "Matrix"), as(y, "Matrix"), p, progress)
  }
)

setMethod(
  "emd_stability",
  signature(
    x = "data.frame",
    y = "data.frame",
    p = "numeric_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(progress)) {
      progress = interactive()
    }
    l = create_matrices_from_data_frames(x, y)
    emd_stability(l$xm, l$ym, p, progress)
  }
)

#' @title Stability of the Embedded Earthmove Distance in Each Sample
#' @description The earthmover distance is the sum of the Minkowski 
#' distances of samples, weighted their transport coefficients. The 
#' `embedded_emd_stability()` function evaluates how stable the distance 
#' is in each of the samples embedding defined by a model by 
#' iteratively removing each sample, calating the new earthmover distance and 
#' subtracting that from the earthmover distance
#' will all samples. The difference between those distances can then be used
#' to evaluate both how dis-simmilar the sample is to other samples as well
#' as how much the distance is dependent on the sample.
#' @aliases
#'   embedded_emd_stability,data.frame,data.frame,rfsrc,numeric_or_missing,logical_or_missing-method 
#' @param x a `matrix`, `Matrix`, or `data.frame`.
#' @param y a `matrix`, `Matrix`, or `data.frame`.
#' @param model the learner that will induce the embedded 
#' @param p an exponent for the order of the distance (default: 2).
#' @param progress Should the progress be shown as the calculation is being
#' performed (default: `interactive()`)? 
#' @return a tibble with the following variable.
#' * var: The variable (either "x" or "y").
#' * row: The row for the corresponding variable.
#' * diff: The difference between the earthmover distance and the earthmover
#'   distance with the row removed.
#' * p_overall: the empirical p-value of the difference compared to all other
#'   differences.
#' * p_within: the empirical p-value of the differences compared to other
#'   difference in the same data set.
#' * p_across: the empirical p-value of the difference compared to differences
#'   in the other data set.
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(randomForestSRC)
#'
#' # Subset iris for our x and y.
#' iris1 = iris |> filter(Sepal.Length < 5.8)
#' iris2 = iris |> filter(Sepal.Length >= 5.8)
#'
#' # Create a model.
#' fit1 = rfsrc(
#'   Species ~ .,
#'   data = iris1
#' )
#'
#' # Calculate the embedding distance between the data sets.
#' embedded_emd_stability(iris1, iris2, fit1)
#' }
#' @docType methods
#' @rdname embedded_emd_stability-methods
#' @export
setGeneric(
  "embedded_emd_stability",
  function(x, y, model, p, progress) standardGeneric("embedded_emd_stability")
)

setMethod(
  "embedded_emd_stability",
  signature(
    x = "data.frame",
    y = "data.frame",
    model = "rfsrc",
    p = "numeric_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, model, p, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(progress)) {
      progress = interactive()
    }
    l = embed_samples(x, y, model)
    emd_stability(l$xm, l$ym, p, progress)
  }
)



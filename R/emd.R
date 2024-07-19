
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
setClassUnion("character_or_missing", c("missing", "character"))
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
#' @param normality_assumtion can we assume the distances are distributed
#' as normal? Possible values are:
#' * `TRUE` assume normality give a warning if the normality test fails.
#' * `FALSE` do not assume normality stop if the normality test fails.
#' Default is `TRUE`.
#' @param progress Should the progress be shown as the calculation is being
#' performed (default: `interactive()`)? 
#' @return return an instance of type `earthmover_stability` holding
#' slots telling whether the distances are assumed normal along with a
#' `tbl_df` with the following variables:
#' * var: The variable (either "x" or "y").
#' * row: The row for the corresponding variable.
#' * dist_across: The earthmover distance from the data set, minus the sample 
#'   to the other data set.
#' * dist_within: The earthmover distance from the data set to itself, minus 
#'   the sample to the other data set.
#' * p_across: the p-value of the sample with respect to the other data set.
#' * p_within: the p-value of the sample with respect to its own data set.
#' @examples
#' X = matrix(rnorm(3*2, mean=-1),ncol=2) # m obs. for X
#' Y = matrix(rnorm(5*2, mean=+1),ncol=2) # n obs. for Y
#' emd_stability(X, Y)
#' @docType methods
#' @rdname emd_stability-methods
#' @export
setGeneric(
  "emd_stability",
  function(x, y, p, normality_assumption, p_value_correction, progress) 
    standardGeneric("emd_stability")
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

left_emd_dist = function(x, y, p, progress) {
  future_map_dbl(
    seq_len(nrow(x)),
    ~ emd_impl(x[-.x,,drop = FALSE], y, p)$dist,
    .options = furrr_options(seed=TRUE),
    .progress = progress
  )
}

setClass(
  "earthmover_stability",
  representation(
    emds = "tbl_df",
    normality = "logical"
  ),
  prototype(
    emds = tibble(),
    normality = NA
  )
)

#' @importFrom furrr future_map_dbl furrr_options
#' @importFrom stats p.adjust
#' @importFrom tibble tibble
setMethod(
  "emd_stability",
  signature(
    x = "Matrix",
    y = "Matrix",
    p = "numeric_or_missing",
    normality_assumption = "logical_or_missing",
    p_value_correction = "character_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, normality_assumption, p_value_correction, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(normality_assumption)) {
      normality_assumption = TRUE
    }
    if (missing(p_value_correction)) {
      p_value_correction = "bonferroni"
    }
    if (missing(progress)) {
      progress = interactive()
    }
    ret = tibble(
      var = c(
        rep("x", nrow(x)),
        rep("y", nrow(y))
      ),
      row = seq_along(var)
    )
    ret$dist_across = c(
      left_emd_dist(x, y, p, progress),
      left_emd_dist(y, x, p, progress)
    )
    ret$dist_within = c(
      left_emd_dist(x, x, p, progress),
      left_emd_dist(y, y, p, progress)
    )
    if (!normality_assumption) {
      # Get the empirical p-values
      perc_fun = \(x, y) {
        r = ecdf(y)(x) 
        vapply(r, \(x) min(x, 1-x), NA_real_) |>
          p.adjust(method = p_value_correction)
      }
    } else {
      xdw = ret$dist_within[ret$var == "x"]
      ksxp = ks.test(xdw, "pnorm", mean(xdw), sd(xdw))
      if (ksxp$p <= 0.05) {
        warning(
          "`x` argument is different from normal with ks p-value of ",
          round(ksxp$p, 3)
        )
      }
      ydw = ret$dist_within[ret$var == "y"]
      ksyp = ks.test(ydw, "pnorm", mean(ydw), sd(ydw))
      if (ksyp$p <= 0.05) {
        warning(
          "`x` argument is different from normal with ks p-value of ",
          round(ksyp$p, 3)
        )
      }
      perc_fun = \(x, y) {
        r = pnorm(x, mean(y), sd(y))
        vapply(r, \(x) min(x, 1-x), NA_real_) |>
          p.adjust(method = p_value_correction)
      }
    }
    ret$p_across = c(
      perc_fun(
        ret$dist_across[ret$var == "x"],
        ret$dist_across[ret$var == "y"]
      ),
      perc_fun(
        ret$dist_across[ret$var == "y"],
        ret$dist_across[ret$var == "x"]
      )
    )
    ret$p_within = c(
      perc_fun(
        ret$dist_across[ret$var == "x"],
        ret$dist_across[ret$var == "x"]
      ),
      perc_fun(
        ret$dist_across[ret$var == "y"],
        ret$dist_across[ret$var == "y"]
      )
    )
    class(ret) = c("earthmover_stability", class(ret))
    new(
      "earthmover_stability",
      emds = tibble(ret),
      normality = normality_assumption
    )
  }
)

#' @importFrom methods as
setMethod(
  "emd_stability",
  signature(
    x = "matrix",
    y = "matrix",
    p = "numeric_or_missing",
    normality_assumption = "logical_or_missing",
    p_value_correction = "character_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, normality_assumption, p_value_correction, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(normality_assumption)) {
      normality_assumption = TRUE
    }
    if (missing(p_value_correction)) {
      p_value_correction = "bonferroni"
    }
    if (missing(progress)) {
      progress = interactive()
    }
    emd_stability(
      as(x, "Matrix"), 
      as(y, "Matrix"), 
      p, 
      normality_assumption, 
      p_value_correction,
      progress
    )
  }
)

setMethod(
  "emd_stability",
  signature(
    x = "data.frame",
    y = "data.frame",
    p = "numeric_or_missing",
    normality_assumption = "logical_or_missing",
    p_value_correction = "character_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, normality_assumption, p_value_correction, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(normality_assumption)) {
      normality_assumption = TRUE
    }
    if (missing(p_value_correction)) {
      p_value_correction = "bonferroni"
    }
    if (missing(progress)) {
      progress = interactive()
    }
    l = create_matrices_from_data_frames(x, y)
    emd_stability(
      l$xm, 
      l$ym, 
      p, 
      normality_assumption, 
      p_value_correction, 
      progress
    )
  }
)



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
  row_outer_product(x, y, \(x, y) sum((x - y)^p)^(1/p))
}

setClassUnion("numeric_or_missing", c("missing", "numeric"))
setClassUnion("logical_or_missing", c("missing", "logical"))
setClassUnion("character_or_missing", c("missing", "character"))

setOldClass(c("earthmover_stability", "ks.test"))

#' @import Matrix
setClassUnion(
  "Matrix_or_matrix_or_data_frame", 
  c("Matrix", "matrix", "data.frame")
)

#' @title The Earthmover Distance
#' @description Calculate the earthmover distance between data sets `x` and 
#' `y`.
#' @aliases 
#'   emd,data.frame,data.frame,numeric_or_missing-method
#'   emd,matrix,matrix,numeric_or_missing-method
#'   emd,Matrix,Matrix,numeric_or_missing-method
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
#' @importFrom tibble as_tibble
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
  ret = list(
    # Log this for numerical stability?
    dist = sum(matches[,3] * apply((xp - yp)^p, 1, sum))^(1/p),
    pairs = matches |> as_tibble()
  )
  class(ret) = "earthmover_dist"
  ret
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
#'   emd_stability,data.frame,data.frame,numeric_or_missing,logical_or_missing,character_or_missing,logical_or_missing-method
#'   emd_stability,matrix,matrix,numeric_or_missing,logical_or_missing,character_or_missing,logical_or_missing-method
#'   emd_stability,Matrix,Matrix,numeric_or_missing,logical_or_missing,character_or_missing,logical_or_missing-method
#' @param x a `matrix`, `Matrix`, or `data.frame`.
#' @param y a `matrix`, `Matrix`, or `data.frame`.
#' @param p an exponent for the order of the distance (default: 2).
#' @param progress Should the progress be shown as the calculation is being
#' performed (default: `FALSE`)? 
#' @return return an instance of type `earthmover_stability` holding
#' slots telling whether the distances are assumed normal along with a
#' `tbl_df` with the following variables:
#' * var: The variable (either "x" or "y").
#' * row: The row for the corresponding variable.
#' * dist_omnibus: The earthmover distance from the data set, minus the sample 
#'   to the other data set.
#' * dist_within: The earthmover distance from the data set to itself, minus 
#'   the sample to the other data set.
#' * p_omnibus: the p-value of the sample with respect to the other data set.
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
  function(x, y, p, adjust, progress) 
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

#' @importFrom tibble tibble
setClass(
  "earthmover_stability",
  slots = c(
    jack_dists = "tbl_df",
    p = "numeric"
  )
)

#' @importFrom dplyr group_by mutate ungroup row_number
#' @importFrom furrr future_map_dbl furrr_options
#' @importFrom methods new
#' @importFrom tibble tibble
#' @importFrom stats p.adjust ks.test pnorm sd 
setMethod(
  "emd_stability",
  signature(
    x = "Matrix",
    y = "Matrix",
    p = "numeric_or_missing",
    adjust = "character_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, adjust, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(adjust)) {
      adjust = "BH"
    }
    if (missing(progress)) {
      progress = FALSE
    }
    ems = tibble(
        var = c(
          rep("x", nrow(x)),
          rep("y", nrow(y))
        )
      ) |>
      group_by(var) |>
      mutate(row = row_number()) |>
      ungroup()

    ems$dist_omnibus = c(
      left_emd_dist(x, y, p, progress),
      left_emd_dist(y, x, p, progress)
    )
    ems$dist_within = c(
      left_emd_dist(x, x, p, progress),
      left_emd_dist(y, y, p, progress)
    )
    dxw = ems$dist_within[ems$var == "x"]^p
    dyw = ems$dist_within[ems$var == "y"]^p
    da = ems$dist_omnibus^p
    ems$p_omnibus = pnorm(da, mean(da), sd(da)) |>
      vapply(\(r) min(r, 1-r), NA_real_) |>
      p.adjust(adjust)
    ems$p_within = c(
      pnorm(dxw, mean(dxw), sd(dxw)) |>
        vapply(\(r) min(r, 1-r), NA_real_) |>
        p.adjust(adjust),
      pnorm(dyw, mean(dyw), sd(dyw)) |>
        vapply(\(r) min(r, 1-r), NA_real_) |>
        p.adjust(adjust)
    )
    new(
      "earthmover_stability",
      jack_dists = ems,
      p = p
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
    adjust = "character_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, adjust, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(adjust)) {
      adjust = "BH"
    }
    if (missing(progress)) {
      progress = FALSE
    }
    emd_stability(
      as(x, "Matrix"), 
      as(y, "Matrix"), 
      p, 
      adjust,
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
    adjust = "character_or_missing",
    progress = "logical_or_missing"
  ),
  function(x, y, p, adjust, progress) {
    if (missing(p)) {
      p = 2
    }
    if (missing(adjust)) {
      adjust = "BH"
    }
    if (missing(progress)) {
      progress = FALSE
    }
    l = create_matrices_from_data_frames(x, y)
    emd_stability(
      l$xm, 
      l$ym, 
      p, 
      adjust,
      progress
    )
  }
)

# @param p_value_correction What type of p-value correction should be 
# employed. Options are c("holm", "hochberg", "hommel", "bonferroni", "BH", 
# "BY", "fdr", "none"). Default is "bonferroni".

# Which points are outliers with respect to their own distribution?
# Which points are outliers with respect to the omnibus distribution.

setClass(
  "earthmover_stability_summary",
  slots = c(
    ems = "earthmover_stability",
    ks = "ks.test",
    within_p_x = "tbl_df",
    within_p_y = "tbl_df",
    omnibus_p = "tbl_df",
    p_val_thresh = "numeric"
  )
)

setMethod(
  "show",
  signature (
    object = "earthmover_stability_summary"
  ),
  function(object) {
    print(object)
  }
)


#' @importFrom cli cli_h1 cli_h3 cli_text style_bold
setMethod(
  "print",
  signature(
    x = "earthmover_stability_summary"
  ),
  function(x, ...) {
    cli_text()
    cli_h1("Data set difference test")
    cli_text()
    cli_text(style_bold("{x@ks$method}"))
    cli_text(
      "{names(x@ks$statistic)} = {round(x@ks$statistic, 3)}, ",
      "p-value = {round(x@ks$p.value, 3)}"
    )
    cli_h1("Within data set outliers")

    cli_h3("Left data set.")

    if (nrow(x@within_p_x) > 0) {
      x@within_p_x |>
        select(row, p_within) |>
        print()
    } else {
      cli_text("(None)")
    }

    cli_h3("Right data set.")

    if (nrow(x@within_p_y) > 0) {
      x@within_p_y |>
        select(row, p_within) |>
        print()
    } else {
      cli_text("(None)")
    }

    cli_h1("Omnibus data set outliers")

    if (nrow(x@omnibus_p) > 0) {
      x@omnibus_p |> 
        mutate(var = if_else(var == "x", "left", "right")) |>
        select(var, row, p_omnibus) |>
        print()
    } else {
      cli_text("(None)")
    }

    cli_text()
  }
)

#' @importFrom dplyr filter select rename if_else
#' @importFrom methods new
setMethod(
  "summary",
  signature(
    object = "earthmover_stability"
  ),
  function(object, p_val_thresh = 0.05, ...) {

    ks = ks.test(
      object@jack_dists$dist_within[object@jack_dists$var == "x"],
      object@jack_dists$dist_within[object@jack_dists$var == "y"]
    )

    within_p_x = object@jack_dists |>
      filter(p_within <= p_val_thresh & var == "x")

    within_p_y = object@jack_dists |>
      filter(p_within <= p_val_thresh & var == "y")


    omnibus_p = object@jack_dists |>
      filter(p_omnibus <= p_val_thresh) 
    new(
      "earthmover_stability_summary",
      ems = object,
      ks = ks,
      within_p_x = within_p_x,
      within_p_y = within_p_y,
      omnibus_p = omnibus_p,
      p_val_thresh = p_val_thresh
    )
  }
)

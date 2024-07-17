
<!-- README.md is generated from README.Rmd. Please edit that file -->

# earthmover

<!-- badges: start -->

[![R-CMD-check](https://github.com/kaneplusplus/earthmover/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kaneplusplus/earthmover/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The earthmover distance calculates distances between sets of data. It is
calculated as the sum of the Minkowski distances of samples, weighted by
their transport coefficient. This package implements the classic
earthmover distance, the earthmover distance for samples that are
embedded by a model, including random forests, as well as analyses and
visualizations to assess the stability of the distances with respect to
their constituent samples.

## Installation

You can install the development version of earthmover from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kaneplusplus/earthmover")
```

## Example

``` r
library(earthmover)
library(future)
library(tibble)

# Run in parallel.
plan("multisession", workers = 2)

# Subset iris for our x and y.
iris1 = iris[1:50,]
iris2 = iris[51:100,]

# Formatting output for nicer presentation.
pairs_to_tibble = \(x) {
  x$pairs = as_tibble(x$pairs, .name_repair = "minimal")
  x
}

# The distance between the two subsets of data.
emd(iris1, iris2) |> 
  pairs_to_tibble()
#> $dist
#> [1] 3.541949
#> 
#> $pairs
#> # A tibble: 50 × 3
#>       ``    ``    ``
#>    <dbl> <dbl> <dbl>
#>  1     1    14  0.02
#>  2     2    19  0.02
#>  3     3    31  0.02
#>  4     4    40  0.02
#>  5     5    47  0.02
#>  6     6    28  0.02
#>  7     7    35  0.02
#>  8     8    18  0.02
#>  9     9    44  0.02
#> 10    10    13  0.02
#> # ℹ 40 more rows

# The distance between the a subsets and itself.
emd(iris1, iris1) |>
  pairs_to_tibble()
#> $dist
#> [1] 0
#> 
#> $pairs
#> # A tibble: 50 × 3
#>       ``    ``    ``
#>    <dbl> <dbl> <dbl>
#>  1     1     1  0.02
#>  2     2     2  0.02
#>  3     3     3  0.02
#>  4     4     4  0.02
#>  5     5     5  0.02
#>  6     6     6  0.02
#>  7     7     7  0.02
#>  8     8     8  0.02
#>  9     9     9  0.02
#> 10    10    10  0.02
#> # ℹ 40 more rows
```

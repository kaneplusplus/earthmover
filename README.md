
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
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(earthmover)
library(future)

# Run in parallel.
plan("multisession", workers = 2)

# Subset iris for our x and y.
iris1 = iris |> filter(Sepal.Length < 5.8)
iris2 = iris |> filter(Sepal.Length >= 5.8)

# The distance between the two subsets of data.

# The distance between the a subsets and itself.
emd(iris1, iris1)
#> $dist
#> [1] 0
#> 
#> $pairs
#>                    
#> 1   1  1 0.01369863
#> 2   2  2 0.01369863
#> 3   3  3 0.01369863
#> 4   4  4 0.01369863
#> 5   5  5 0.01369863
#> 6   6  6 0.01369863
#> 7   7  7 0.01369863
#> 8   8  8 0.01369863
#> 9   9  9 0.01369863
#> 10 10 10 0.01369863
#> 11 11 11 0.01369863
#> 12 12 12 0.01369863
#> 13 13 13 0.01369863
#> 14 14 14 0.01369863
#> 15 15 15 0.01369863
#> 16 16 16 0.01369863
#> 17 17 17 0.01369863
#> 18 18 18 0.01369863
#> 19 19 19 0.01369863
#> 20 20 20 0.01369863
#> 21 21 21 0.01369863
#> 22 22 22 0.01369863
#> 23 23 23 0.01369863
#> 24 24 24 0.01369863
#> 25 25 25 0.01369863
#> 26 26 26 0.01369863
#> 27 27 27 0.01369863
#> 28 28 28 0.01369863
#> 29 29 29 0.01369863
#> 30 30 30 0.01369863
#> 31 31 31 0.01369863
#> 32 32 32 0.01369863
#> 33 33 33 0.01369863
#> 34 34 34 0.01369863
#> 35 35 35 0.01369863
#> 36 36 36 0.01369863
#> 37 37 37 0.01369863
#> 38 38 38 0.01369863
#> 39 39 39 0.01369863
#> 40 40 40 0.01369863
#> 41 41 41 0.01369863
#> 42 42 42 0.01369863
#> 43 43 43 0.01369863
#> 44 44 44 0.01369863
#> 45 45 45 0.01369863
#> 46 46 46 0.01369863
#> 47 47 47 0.01369863
#> 48 48 48 0.01369863
#> 49 49 49 0.01369863
#> 50 50 50 0.01369863
#> 51 51 51 0.01369863
#> 52 52 52 0.01369863
#> 53 53 53 0.01369863
#> 54 54 54 0.01369863
#> 55 55 55 0.01369863
#> 56 56 56 0.01369863
#> 57 57 57 0.01369863
#> 58 58 58 0.01369863
#> 59 59 59 0.01369863
#> 60 60 60 0.01369863
#> 61 61 61 0.01369863
#> 62 62 62 0.01369863
#> 63 63 63 0.01369863
#> 64 64 64 0.01369863
#> 65 65 65 0.01369863
#> 66 66 66 0.01369863
#> 67 67 67 0.01369863
#> 68 68 68 0.01369863
#> 69 69 69 0.01369863
#> 70 70 70 0.01369863
#> 71 71 71 0.01369863
#> 72 72 72 0.01369863
#> 73 73 73 0.01369863
```

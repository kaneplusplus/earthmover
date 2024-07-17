---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# earthmover

<!-- badges: start -->
[![R-CMD-check](https://github.com/kaneplusplus/earthmover/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kaneplusplus/earthmover/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The earthmover distance calculates distances between sets of data. It is 
calculated as the sum of the Minkowski distances of samples, weighted by 
their transport coefficient. This package implements the classic earthmover 
distance, the earthmover distance for samples that are embedded by a model, 
including random forests, as well as analyses and visualizations to assess the 
stability of the distances with respect to their constituent samples.

## Installation

You can install the development version of earthmover from [GitHub](https://github.com/) with:

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
#> Error in library(earthmover): there is no package called 'earthmover'
library(future)

# Run in parallel.
plan("multisession", workers = 2)

# Subset iris for our x and y.
iris1 = iris |> filter(Sepal.Length < 5.8)
iris2 = iris |> filter(Sepal.Length >= 5.8)

# The distance between the two subsets of data.

# The distance between the a subsets and itself.
emd(iris1, iris1)
#> Error in emd(iris1, iris1): could not find function "emd"
```


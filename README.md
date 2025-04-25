
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package NCor

<!-- badges: start -->
<!-- badges: end -->

## Overview

This package accompanies the paper “Proper Correlation Coefficients for
Nominal Random Variables” by Jan-Lukas Wermuth.

It provides an estimator for the coefficient proposed in the paper
together with a confidence interval and a p-value for the corresponding
independence test.

## Installation

You can install the development version of NCor from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
devtools::install_github("jan-lukas-wermuth/NCor")
devtools::install_github("jan-lukas-wermuth/RCor")
```

## Examples

The first example simulates a bivariate sample with one nominal and one
continuous marginal distribution that are independent if and only if
`rho = 0`.

``` r
library(NCor)
library(RCor)

# Define parameters
n <- 100 # sample size
rho <- 1 # dependence parameter

# Generate nominal covariate and create continuous dependent variable via regression
set.seed(1) 
X <- sample(c("A", "B", "C"), n, prob = c(1/3, 1/3, 1/3), replace = TRUE)
C_ind <- ifelse(X == "C", 1, 0)
B_ind <- ifelse(X == "B", 1, 0)
Y <- rho * B_ind - rho * C_ind + rnorm(n)

# Compute coefficient together with an estimated 90% confidence interval and the p-value for the related independence test
NCor::NCor(X, Y, nominal = "r", CIs = TRUE, Test = TRUE)[[1]]
#> # A tibble: 1 × 4
#>   Gamma CI_lower CI_upper   PValue
#>   <dbl>    <dbl>    <dbl>    <dbl>
#> 1 0.581    0.471    0.692 7.64e-10
```

A second example simulates from a $3\times 3$ probability contingency
table that has one skewed and one uniform marginal distribution which
are also independent if and only if `rho = 0`.

``` r
library(NCor)
library(RCor)

# Define a 3x3 probability contingency table with dependence
rho <- 0.01 # dependence parameter
probabilities <- matrix(c(76/300 + 2*rho, 76/300 + 2*rho, 76/300 - 4*rho, 4/100 - rho, 4/100 - rho, 4/100 + 2*rho, 4/100 - rho, 4/100 - rho, 4/100 + 2*rho), nrow = 3)
prob_vector <- as.vector(probabilities)

# Simulate data from the probability contingency table
n <- 100
set.seed(1)
sampled_indices <- sample(1:9, size = n, replace = TRUE, prob = prob_vector)
contingency_table <- table(factor(sampled_indices, levels = 1:9))
contingency_matrix <- matrix(contingency_table, nrow = 3)

# Compute coefficient together with an estimated 90% confidence interval and the p-value for the related independence test
NCor::NCor(contingency_matrix, nominal = "rc", CIs = TRUE, Test = TRUE)[[1]]
#> # A tibble: 1 × 4
#>   Gamma CI_lower CI_upper  PValue
#>   <dbl>    <dbl>    <dbl>   <dbl>
#> 1 0.611    0.395    0.828 0.00263
```

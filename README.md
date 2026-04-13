# mqriskR

<!-- badges: start -->
<!-- badges: end -->

The goal of mqriskR is to provide functions for actuarial risk modeling,
including survival models, life annuities, multiple-decrement models,
and mortality improvement projections, using standard actuarial notation.

The package is designed to align with standard actuarial notation and
supports teaching, exam preparation, and reproducible actuarial analysis.

## Installation

You can install the development version of mqriskR from GitHub:

```r
# install.packages("remotes")
remotes::install_github("aokine/mqriskR")
```

## Example

This example computes a whole life annuity under a uniform distribution of deaths.

``` r
library(mqriskR)

ax(40, i = 0.05, model = "uniform", omega = 100)
```


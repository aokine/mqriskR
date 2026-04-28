# mqriskR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mqriskR)](https://CRAN.R-project.org/package=mqriskR)
<!-- badges: end -->

`mqriskR` provides functions for actuarial risk modeling, including survival models, life annuities, multiple-decrement models, and mortality improvement projections, using standard actuarial notation.

The package is designed to support teaching, exam preparation, and reproducible actuarial analysis.

## Installation

You can install the stable version from CRAN:

```r
install.packages("mqriskR")
```

You can install the development version from GitHub:

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

## Features

- Life insurance and annuity functions (discrete, continuous, and m-thly)
- Survival models and mortality laws
- Multiple-decrement models
- Mortality improvement projections
- Functions aligned with standard actuarial notation

## Purpose

This package is intended for:
	
- Actuarial students preparing for professional exams
- Instructors teaching life contingencies
- Practitioners needing transparent and reproducible calculations


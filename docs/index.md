# mqriskR

`mqriskR` provides functions for actuarial risk modeling, including
survival models, life annuities, multiple-decrement models, and
mortality improvement projections, using standard actuarial notation.

The package is designed to support teaching, exam preparation, and
reproducible actuarial analysis.

## Installation

You can install the stable version from CRAN:

``` r

install.packages("mqriskR")
```

You can install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("aokine/mqriskR")
```

## Example

This example computes a whole life annuity under a uniform distribution
of deaths.

``` r

library(mqriskR)

ax(40, i = 0.05, model = "uniform", omega = 100)
```

## Features

- Life insurance and annuity functions (discrete, continuous, and
  m-thly)
- Survival models and mortality laws
- Multiple-decrement models
- Mortality improvement projections
- Functions aligned with standard actuarial notation

## Purpose

This package is intended for:

- Actuarial students preparing for professional exams
- Instructors teaching life contingencies
- Practitioners needing transparent and reproducible calculations

## Documentation

Full function documentation is available on the pkgdown site:

<https://aokine.github.io/mqriskR/>

## References

The methods implemented in this package are aligned with standard
actuarial texts, including:

- Camilli, S. J., Duncan, I., and London, R. L. (2014,
  <ISBN:9781625423474>) “Models for Quantifying Risk”, 6th Edition,
  ACTEX Publications.

- Dickson, D. C. M., Hardy, M. R., and Waters, H. R. (2020,
  <ISBN:9781108478083>) “Actuarial Mathematics for Life Contingent
  Risks”.

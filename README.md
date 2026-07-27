# mqriskR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mqriskR)](https://CRAN.R-project.org/package=mqriskR)
<!-- badges: end -->

**mqriskR** is an R package for actuarial mathematics and life contingency modeling. It provides functions for calculating actuarial present values, premiums, reserves, survival probabilities, life annuities, pension benefits, multiple-decrement models, and mortality improvement using both life tables and parametric survival models.

The package is designed to support actuarial education, professional exam preparation, research, and reproducible actuarial analysis.

## Installation

Install the stable release from CRAN:

```r
install.packages("mqriskR")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("aokine/mqriskR")
```

## Example

Compute the actuarial present value of a whole life annuity under a uniform distribution of deaths.

```r
library(mqriskR)

ax(
  x = 40,
  i = 0.05,
  model = "uniform",
  omega = 100
)
```

## Main Features

- Life insurance present values
- Life annuities (discrete, continuous, and m-thly)
- Premium calculations
- Policy reserve calculations
- Pension mathematics
- Survival probabilities
- Multiple-life models
- Multiple-decrement models
- Mortality improvement projections
- Spot interest rate models
- Variable interest models
- Support for both life tables and parametric survival models
- Functions using standard actuarial notation

## Version 0.1.1

Version 0.1.1 focuses on improving the quality, consistency, and usability of the package.

Highlights include:

- Improved support for life table objects
- Expanded vectorized input support
- Improved handling of finite life tables
- More consistent behavior across related functions
- Cleaner, descriptive documentation
- Improved examples and package manual
- Enhanced input validation and error handling

No breaking changes were introduced. Existing code written for earlier versions of **mqriskR** continues to work.

## Intended Audience

The package is intended for:

- Actuarial students preparing for professional examinations
- University instructors teaching actuarial mathematics
- Researchers developing actuarial methods
- Practicing actuaries requiring transparent and reproducible calculations

## Documentation

Complete documentation for all exported functions is available on the pkgdown website:

https://aokine.github.io/mqriskR/

## References

The methods implemented in **mqriskR** are based on standard actuarial references, including:

- Camilli, S. J., Duncan, I., and London, R. L. (2014). *Models for Quantifying Risk* (6th ed.). ACTEX Publications.

- Dickson, D. C. M., Hardy, M. R., and Waters, H. R. (2020). *Actuarial Mathematics for Life Contingent Risks* (2nd ed.). Cambridge University Press.

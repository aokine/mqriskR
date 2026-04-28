# Variance of whole life insurance PV

Computes \\\mathrm{Var}(Z_x) = {}^{2}A_x - A_x^2\\.

## Usage

``` r
var_Ax(x, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to survival-model functions.

- tol:

  Numerical tolerance for truncating infinite sums.

- k_max:

  Maximum number of terms in the sum.

## Value

Numeric vector of variances.

# Variance of deferred insurance PV

Variance of deferred insurance PV

## Usage

``` r
var_nAx(x, n, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000)
```

## Arguments

- x:

  Age.

- n:

  Deferral period.

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

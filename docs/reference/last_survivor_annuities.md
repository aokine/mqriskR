# Last-survivor annuity functions

Computes temporary and whole-life annuities for two lives under the
last-survivor status.

## Usage

``` r
adotxybarn(x, y, n, i, tbl = NULL, model = NULL, ...)

axybarn(x, y, n, i, tbl = NULL, model = NULL, ...)

adotxybar(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)

axybar(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

- n:

  Term in years. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the survival model or life-table
  functions.

- k_max:

  Maximum number of terms used for an infinite series.

- tol:

  Numerical tolerance used to assess convergence.

## Value

A numeric vector of actuarial present values.

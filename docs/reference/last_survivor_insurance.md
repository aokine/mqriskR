# Last-survivor insurance functions

Computes temporary and whole-life insurance values for two lives under
the last-survivor status.

## Usage

``` r
Axybarn1(x, y, n, i, tbl = NULL, model = NULL, ...)

Axybarn(x, y, n, i, tbl = NULL, model = NULL, ...)

Axybar(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)
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

## Details

`Axybarn1()` computes temporary insurance payable at the second death
within the term.

`Axybarn()` computes last-survivor endowment insurance.

`Axybar()` computes last-survivor whole-life insurance.

# Joint-life insurance functions

Computes temporary and whole-life insurance values for two lives under
the joint-life status.

## Usage

``` r
Axyn1(x, y, n, i, tbl = NULL, model = NULL, ...)

Axyn(x, y, n, i, tbl = NULL, model = NULL, ...)

Axy(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)
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

`Axyn1()` computes an `n`-year joint-life term insurance.

`Axyn()` computes an `n`-year joint-life endowment insurance.

`Axy()` computes joint-life whole-life insurance.

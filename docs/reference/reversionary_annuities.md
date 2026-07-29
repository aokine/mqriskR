# Reversionary annuity functions

Computes reversionary whole-life annuities payable to one life after the
death of the other life.

## Usage

``` r
ax_y(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)

ay_x(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

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

`ax_y()` computes an annuity payable to the second life after the death
of the first life.

`ay_x()` computes an annuity payable to the first life after the death
of the second life.

# Reversionary annuity to (x) after death of (y)

Computes \\a\_{y\|x} = a_x - a\_{xy}\\.

## Usage

``` r
ay_x(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12)
```

## Arguments

- x:

  Age of first life.

- y:

  Age of second life.

- i:

  Effective annual interest rate.

- tbl:

  Life table.

- model:

  Survival model.

- ...:

  Additional model parameters.

- k_max:

  Maximum number of terms.

- tol:

  Convergence tolerance.

## Value

Numeric vector.

## Examples

``` r
ay_x(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 3.413213
```

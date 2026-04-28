# Joint-life whole life annuity-due

Computes \\\ddot{a}\_{xy}\\.

## Usage

``` r
adotxy(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12)
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
adotxy(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 10.96154
```

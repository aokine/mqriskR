# Last-survivor whole life insurance

Computes \\A\_{\overline{xy}} = A_x + A_y - A\_{xy}\\.

## Usage

``` r
Axybar(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12)
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
Axybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.2025846
```

# Joint-life temporary annuity-due

Computes \\\ddot{a}\_{xy:\overline{n}\|}\\.

## Usage

``` r
adotxyn(x, y, n, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of first life.

- y:

  Age of second life.

- n:

  Term.

- i:

  Effective annual interest rate.

- tbl:

  Life table.

- model:

  Survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

## Examples

``` r
adotxyn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 6.956659
```

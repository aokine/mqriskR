# Joint-life temporary annuity-immediate

Computes \\a\_{xy:\overline{n}\|}\\.

## Usage

``` r
axyn(x, y, n, i, tbl = NULL, model = NULL, ...)
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

## Details

Shared documentation topic used to avoid filename collisions with
case-distinct function names on case-insensitive file systems.

## Examples

``` r
axyn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 6.547384
```

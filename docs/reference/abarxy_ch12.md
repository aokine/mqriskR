# Continuous joint-life whole life annuity

Computes \\\overline{a}\_{xy} = \int_0^\infty v^t\\{}\_tp\_{xy}\\dt\\.

## Usage

``` r
abarxy(x, y, i, tbl = NULL, model = NULL, ...)
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

## Value

Numeric vector.

## Details

Shared documentation topic used to avoid filename collisions with
case-distinct function names on case-insensitive file systems.

## Examples

``` r
abarxy(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 10.45444
```

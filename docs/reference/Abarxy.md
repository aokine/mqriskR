# Continuous joint-life whole life insurance

Computes \\\overline{A}\_{xy} = 1 - \delta \overline{a}\_{xy}\\.

## Usage

``` r
Abarxy(x, y, i, tbl = NULL, model = NULL, ...)
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

## Examples

``` r
Abarxy(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.4899262
```

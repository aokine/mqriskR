# Continuous reversionary annuity to (x) after death of (y)

Computes \\\overline{a}\_{y\|x} = \overline{a}\_x -
\overline{a}\_{xy}\\.

## Usage

``` r
abary_x(x, y, i, tbl = NULL, model = NULL, ...)
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
abary_x(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 3.41493
```

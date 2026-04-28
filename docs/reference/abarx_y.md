# Continuous reversionary annuity to (y) after death of (x)

Computes \\\overline{a}\_{x\|y} = \overline{a}\_y -
\overline{a}\_{xy}\\.

## Usage

``` r
abarx_y(x, y, i, tbl = NULL, model = NULL, ...)
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
abarx_y(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 2.372486
```

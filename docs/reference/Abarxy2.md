# Continuous contingent insurance: benefit on death of (x) if after (y)

Computes \\\overline{A}\_{xy}^{2} = \overline{A}\_x -
\overline{A}\_{xy}^{1}\\.

## Usage

``` r
Abarxy2(x, y, i, tbl = NULL, model = NULL, ...)
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
Abarxy2(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.1095292
```

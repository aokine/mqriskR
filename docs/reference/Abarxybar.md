# Continuous last-survivor whole life insurance

Computes \\\overline{A}\_{\overline{xy}} = \overline{A}\_x +
\overline{A}\_y - \overline{A}\_{xy}\\.

## Usage

``` r
Abarxybar(x, y, i, tbl = NULL, model = NULL, ...)
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
Abarxybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.2075573
```

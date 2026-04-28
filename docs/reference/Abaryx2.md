# Continuous contingent insurance: benefit on death of (y) if after (x)

Computes \\\overline{A}\_{xy}^{\hspace{1mm}2} = \overline{A}\_y -
\overline{A}\_{xy}^{\hspace{1mm}1}\\.

## Usage

``` r
Abaryx2(x, y, i, tbl = NULL, model = NULL, ...)
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
Abaryx2(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.09802813
```

# Continuous contingent insurance: benefit on death of (y) if before (x)

Computes \\\overline{A}\_{xy}^{\hspace{1mm}1} = \int_0^\infty
v^t\\{}\_tp\_{xy}\mu\_{y+t}\\dt\\.

## Usage

``` r
Abaryx1(x, y, i, tbl = NULL, model = NULL, ...)
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
Abaryx1(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.2761441
```

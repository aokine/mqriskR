# Probability that (x) fails before (y) within n years

Computes \\{}\_n q\_{xy}^{1} = \int_0^n {}\_tp\_{xy}\mu\_{x+t}\\dt\\
under independence.

## Usage

``` r
tqxy1(x, y, n, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of first life.

- y:

  Age of second life.

- n:

  Term in years.

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
tqxy1(40, 50, n = 10, model = "uniform", omega = 100)
#> [1] 0.15
```

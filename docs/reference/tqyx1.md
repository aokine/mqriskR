# Probability that (y) fails before (x) within n years

Computes \\{}\_n q\_{xy}^{\hspace{1mm}1} = \int_0^n
{}\_tp\_{xy}\mu\_{y+t}\\dt\\ under independence.

## Usage

``` r
tqyx1(x, y, n, tbl = NULL, model = NULL, ...)
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
tqyx1(40, 50, n = 10, model = "uniform", omega = 100)
#> [1] 0.1833333
```

# Probability that (x) fails after (y) within n years

Computes \\{}\_n q\_{xy}^{2} = {}\_n q_x - {}\_n q\_{xy}^{1}\\.

## Usage

``` r
tqxy2(x, y, n, tbl = NULL, model = NULL, ...)
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
tqxy2(40, 50, n = 10, model = "uniform", omega = 100)
#> [1] 0.01666667
```

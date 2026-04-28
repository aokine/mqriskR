# Joint-life failure probability

Computes \\{}\_tq\_{xy} = 1 - {}\_tp\_{xy}\\.

## Usage

``` r
tqxy(x, y, t, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of first life.

- y:

  Age of second life.

- t:

  Time.

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
tqxy(40, 50, t = 10, model = "uniform", omega = 100)
#> [1] 0.3333333
```

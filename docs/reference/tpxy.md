# Joint-life survival probability

Computes \\{}\_tp\_{xy} = {}\_tp_x\\{}\_tp_y\\ under independence.

## Usage

``` r
tpxy(x, y, t, tbl = NULL, model = NULL, ...)
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
tpxy(40, 50, t = 10, model = "uniform", omega = 100)
#> [1] 0.6666667
```

# Last-survivor survival probability

Computes \\{}\_tp\_{\overline{xy}} = {}\_tp_x + {}\_tp_y -
{}\_tp\_{xy}\\.

## Usage

``` r
tpxybar(x, y, t, tbl = NULL, model = NULL, ...)
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
tpxybar(40, 50, t = 10, model = "uniform", omega = 100)
#> [1] 0.9666667
```

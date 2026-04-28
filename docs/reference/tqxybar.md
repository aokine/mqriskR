# Last-survivor failure probability

Computes \\{}\_tq\_{\overline{xy}} = 1 - {}\_tp\_{\overline{xy}}\\.

## Usage

``` r
tqxybar(x, y, t, tbl = NULL, model = NULL, ...)
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
tqxybar(40, 50, t = 10, model = "uniform", omega = 100)
#> [1] 0.03333333
```

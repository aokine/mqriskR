# Joint-life pure endowment

Computes \\{}\_nE\_{xy} = v^n\\{}\_np\_{xy}\\.

## Usage

``` r
nExy(x, y, n, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of first life.

- y:

  Age of second life.

- n:

  Term.

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
nExy(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.4092755
```

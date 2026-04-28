# Fractional-duration endowment reserve

Computes the interpolated reserve \\{}\_{t+s}V = ({}\_tV + P)(1-s) +
{}\_{t+1}V \cdot s\\ for an n-year endowment insurance.

Computes the interpolated reserve \\{}\_{t+s}V = ({}\_tV + P)(1-s) +
{}\_{t+1}V \cdot s\\ for an n-year endowment insurance.

## Usage

``` r
tsVxn(x, n, t, s, i, tbl = NULL, model = NULL, ...)

tsVxn(x, n, t, s, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age.

- n:

  Term in years.

- t:

  Integer duration with t \< n.

- s:

  Fractional part in \[0, 1\].

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

Numeric vector.

## Examples

``` r
tsVxn(40, n = 20, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.388849
tsVxn(40, n = 20, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.388849
```

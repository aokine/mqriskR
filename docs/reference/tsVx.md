# Fractional-duration whole life reserve

Computes the interpolated reserve \\{}\_{t+s}V = ({}\_tV + P_x)(1-s) +
{}\_{t+1}V \cdot s\\ for \\0 \le s \le 1\\.

Computes the interpolated reserve \\{}\_{t+s}V = ({}\_tV + P_x)(1-s) +
{}\_{t+1}V \cdot s\\ for \\0 \le s \le 1\\.

## Usage

``` r
tsVx(x, t, s, i, tbl = NULL, model = NULL, ...)

tsVx(x, t, s, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age.

- t:

  Integer contract duration.

- s:

  Fractional part of duration in \[0, 1\].

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
tsVx(40, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.08762133
tsVx(40, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.08762133
```

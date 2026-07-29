# Whole life reserve with m-thly premiums

Computes the reserve for a whole life insurance funded by true m-thly
premiums.

## Usage

``` r
tVx_m(x, t, m, i, model = NULL, ..., tbl = NULL)
```

## Arguments

- x:

  Issue age.

- t:

  Duration.

- m:

  Number of premium payments per year.

- i:

  Effective annual interest rate.

- model:

  Optional parametric survival model name.

- ...:

  Additional model parameters.

- tbl:

  Optional life table object.

## Value

A numeric vector of values.

## Examples

``` r
tVx_m(40, t = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.07327183
```

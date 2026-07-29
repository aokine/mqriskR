# Whole life reserve with continuous premiums

Computes the reserve for a discrete whole life insurance funded by
continuous premiums.

## Usage

``` r
tVbarx(x, t, i, model = NULL, ..., tbl = NULL)
```

## Arguments

- x:

  Issue age.

- t:

  Duration.

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
tVbarx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.07334295
```

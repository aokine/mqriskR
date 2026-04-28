# Whole life reserve with continuous premiums

Computes the Chapter 10 reserve for a discrete whole life insurance
funded by continuous premiums.

## Usage

``` r
tVbarx(x, t, i, model, ...)
```

## Arguments

- x:

  Issue age.

- t:

  Duration.

- i:

  Effective annual interest rate.

- model:

  Survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

## Examples

``` r
tVbarx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.07334295
```

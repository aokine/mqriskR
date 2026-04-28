# Whole life reserve with m-thly premiums

Computes the Chapter 10 reserve for a whole life insurance funded by
true m-thly premiums.

## Usage

``` r
tVx_m(x, t, m, i, model, ...)
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

  Survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

## Examples

``` r
tVx_m(40, t = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.07327183
```

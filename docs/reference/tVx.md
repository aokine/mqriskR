# Whole life net level premium reserve

Computes the Chapter 10 prospective reserve for a whole life insurance
with annual premiums: reserve at duration t equals future APV of
benefits minus future APV of net premiums.

## Usage

``` r
tVx(x, t, i, model, ...)
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
tVx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.07250474
```

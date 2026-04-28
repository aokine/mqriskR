# h-pay whole life net level premium reserve

Computes the Chapter 10 prospective reserve for an h-pay whole life
policy.

## Usage

``` r
htVx(x, h, t, i, model, ...)
```

## Arguments

- x:

  Issue age.

- h:

  Premium-paying period in years.

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
htVx(40, h = 10, t = 5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.1554969
```

# Deferred annuity-immediate reserve

Computes the reserve for an n-year deferred annuity-immediate for \\t \<
n\\.

## Usage

``` r
tVnax(x, n, t, i, model, ...)
```

## Arguments

- x:

  Issue age.

- n:

  Deferral period.

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
tVnax(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 3.589045
```

# Deferred annuity-immediate premium

Computes \\P({}\_{n\|}a_x) = {}\_{n\|}a_x /
\ddot{a}\_{x:\overline{n}\|}\\.

## Usage

``` r
Pnax(x, n, i, model, ...)
```

## Arguments

- x:

  Issue age.

- n:

  Deferral period.

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
Pnax(40, n = 20, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.2430708
```

# Term insurance net level premium reserve

Computes the Chapter 10 prospective reserve for an n-year term
insurance.

## Usage

``` r
tVxn1(x, n, t, i, model, ...)
```

## Arguments

- x:

  Issue age.

- n:

  Term in years.

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
tVxn1(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.01836757
```

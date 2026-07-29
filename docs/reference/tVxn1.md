# Term insurance net level premium reserve

Computes the prospective reserve for an n-year term insurance.

## Usage

``` r
tVxn1(x, n, t, i, model = NULL, ..., tbl = NULL)
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

  Optional parametric survival model name.

- ...:

  Additional model parameters.

- tbl:

  Optional life table object.

## Value

A numeric vector of values.

## Examples

``` r
tVxn1(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.01836757
```

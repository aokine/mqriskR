# Term insurance reserve by retrospective method

Computes the retrospective term reserve for \\t \le n\\.

## Usage

``` r
tVxn1_ret(x, n, t, i, model, ...)
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
tVxn1_ret(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.01836757
```

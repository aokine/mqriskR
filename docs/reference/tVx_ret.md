# Whole life net level premium reserve by retrospective method

Computes the Chapter 10 retrospective reserve \\{}\_tV_x = P_x
\ddot{s}\_{x:\overline{t}\|} - {}\_tk_x\\.

## Usage

``` r
tVx_ret(x, t, i, model, ...)
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
tVx_ret(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.07250474
```

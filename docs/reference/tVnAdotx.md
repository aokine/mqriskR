# Deferred annuity-due reserve

Computes the Chapter 10 reserve \\{}\_tV({}\_{n\|}\ddot{a}\_x)\\ for \\t
\< n\\.

## Usage

``` r
tVnAdotx(x, n, t, i, model, ...)
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
tVnAdotx(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 3.915575
```

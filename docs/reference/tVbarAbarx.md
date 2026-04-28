# Fully continuous whole life reserve

Computes the Chapter 10 reserve for a whole life insurance with
continuous premiums and immediate payment of claims.

## Usage

``` r
tVbarAbarx(x, t, i, model, ...)
```

## Arguments

- x:

  Issue age.

- t:

  Duration, allowed to be any nonnegative numeric value.

- i:

  Effective annual interest rate.

- model:

  Survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

## Details

In the fully continuous setting, reserve time \\t\\ may be any
nonnegative real value.

## Examples

``` r
tVbarAbarx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.07516161
tVbarAbarx(40, t = c(19, 19.25, 19.5, 19.75, 20), i = 0.06,
           model = "uniform", omega = 100)
#> [1] 0.1422970 0.1447551 0.1472088 0.1496845 0.1521824
```

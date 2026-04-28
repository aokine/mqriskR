# Mean reserve for whole life insurance

Computes the mean reserve \\{}\_{t+1/2}V\\.

Computes the mean reserve \\{}\_{t+1/2}V\\.

## Usage

``` r
meanVx(x, t, i, tbl = NULL, model = NULL, ...)

meanVx(x, t, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age.

- t:

  Integer contract duration.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

Numeric vector.

## Examples

``` r
meanVx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.08762133
meanVx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.08762133
```

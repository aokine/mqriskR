# Full preliminary term first-year modified premium

Computes \\\alpha^F = v q_x = A\_{x:\overline{1}\|}^{1}\\.

Computes \\\alpha^F = v q_x = A\_{x:\overline{1}\|}^{1}\\.

## Usage

``` r
alphaF(x, i, tbl = NULL, model = NULL, ...)

alphaF(x, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age.

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
alphaF(40, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.01587302
alphaF(40, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.01587302
```

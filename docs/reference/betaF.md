# Full preliminary term renewal modified premium

Computes the FPT renewal premium \\\beta^F = P\_{x+1}\\ for whole life
insurance.

Computes the FPT renewal premium \\\beta^F = P\_{x+1}\\ for whole life
insurance.

## Usage

``` r
betaF(x, i, tbl = NULL, model = NULL, ...)

betaF(x, i, tbl = NULL, model = NULL, ...)
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
betaF(40, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.02240155
betaF(40, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.02240155
```

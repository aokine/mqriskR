# Variance of continuous whole life insurance PV

Computes \\\mathrm{Var}(\bar{Z}\_x) = {}^{2}\bar{A}\_x - \bar{A}\_x^2\\.

## Usage

``` r
var_Abarx(x, i, model, ...)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

## Value

Numeric vector of variances.

# Variance of term insurance PV

Computes \\\mathrm{Var}(Z\_{x:\overline{n}\|}^{1}) =
{}^{2}A\_{x:\overline{n}\|}^{1} - (A\_{x:\overline{n}\|}^{1})^2\\.

## Usage

``` r
var_Axn1(x, n, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age.

- n:

  Term.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to survival-model functions.

## Value

Numeric vector of variances.

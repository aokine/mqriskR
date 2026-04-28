# Variance of continuous term insurance PV

Computes \\\mathrm{Var}(\bar{Z}\_{x:\overline{n}\|}^{1}) =
{}^{2}\bar{A}\_{x:\overline{n}\|}^{1} -
(\bar{A}\_{x:\overline{n}\|}^{1})^2\\.

## Usage

``` r
var_Abarxn1(x, n, i, model, ...)
```

## Arguments

- x:

  Age.

- n:

  Term.

- i:

  Effective annual interest rate.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

## Value

Numeric vector of variances.

# Continuous term insurance APV

Computes \\\bar{A}\_{x:\overline{n}\|}^{1} = \int_0^n v^t \\ {}\_t p_x
\mu\_{x+t}\\dt\\.

## Usage

``` r
Abarxn1(x, n, i, model, ...)
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

Numeric vector of APVs.

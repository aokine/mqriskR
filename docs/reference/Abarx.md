# Continuous whole life insurance APV

Computes \\\bar{A}\_x = \int_0^\infty v^t \\ {}\_t p_x \mu\_{x+t}\\dt\\.

## Usage

``` r
Abarx(x, i, model, ...)
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

Numeric vector of APVs.

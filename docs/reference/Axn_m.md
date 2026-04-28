# m-thly endowment insurance APV

Computes \\A\_{x:\overline{n}\|}^{(m)} = A\_{x:\overline{n}\|}^{1(m)} +
v^n\\{}\_np_x\\.

## Usage

``` r
Axn_m(x, n, i, m, model, ...)
```

## Arguments

- x:

  Age.

- n:

  Term.

- i:

  Effective annual interest rate.

- m:

  Positive integer payment frequency.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

## Value

Numeric vector of APVs.

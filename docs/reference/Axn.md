# Endowment insurance APV

Computes \\A\_{x:\overline{n}\|} = A\_{x:\overline{n}\|}^{1} +
{}\_nE_x\\.

## Usage

``` r
Axn(x, n, i, tbl = NULL, model = NULL, ...)
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

Numeric vector of APVs.

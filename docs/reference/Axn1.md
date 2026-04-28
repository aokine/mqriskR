# Term insurance APV

Computes \\A\_{x:\overline{n}\|}^{1} = \sum\_{k=0}^{n-1} v^{k+1}
{}\_{k\mid}q_x\\.

## Usage

``` r
Axn1(x, n, i, tbl = NULL, model = NULL, ...)
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

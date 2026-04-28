# Second moment of pure endowment PV

Computes \\{}^{2}{}\_nE_x = (v')^n {}\_n p_x\\.

## Usage

``` r
A2nEx(x, n, i, tbl = NULL, model = NULL, ...)
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

Numeric vector of second moments.

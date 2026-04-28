# Whole life insurance APV

Computes \\A_x = \sum\_{k=0}^\infty v^{k+1} {}\_{k\mid}q_x\\.

## Usage

``` r
Ax(x, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to survival-model functions.

- tol:

  Numerical tolerance for truncating infinite sums.

- k_max:

  Maximum number of terms in the sum.

## Value

Numeric vector of APVs.

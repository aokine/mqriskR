# m-thly whole life insurance APV

Computes \\A_x^{(m)} = \sum\_{j=0}^\infty v^{(j+1)/m}\Pr(j/m \< T_x \le
(j+1)/m)\\.

## Usage

``` r
Ax_m(x, i, m, model, ..., tol = 1e-12, j_max = 100000L)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- m:

  Positive integer payment frequency.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

- tol:

  Numerical tolerance for truncating the infinite sum.

- j_max:

  Maximum number of m-thly intervals in the sum.

## Value

Numeric vector of APVs.

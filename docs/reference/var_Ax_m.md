# Variance of m-thly whole life insurance PV

Computes \\\mathrm{Var}(Z_x^{(m)}) = {}^{2}A_x^{(m)} - (A_x^{(m)})^2\\.

## Usage

``` r
var_Ax_m(x, i, m, model, ..., tol = 1e-12, j_max = 100000L)
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

Numeric vector of variances.

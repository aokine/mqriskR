# Increasing whole life insurance

Computes \$\$(IA)\_x = \sum\_{k=0}^{\infty} (k+1) v^{k+1} \Pr(K_x =
k).\$\$

## Usage

``` r
IAx(x, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional model parameters.

- tol:

  Numerical tolerance for truncation.

- k_max:

  Maximum number of terms.

## Value

Numeric vector.

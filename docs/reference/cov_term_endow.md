# Covariance of term insurance and pure endowment PVs

Computes \\\mathrm{Cov}(Z\_{x:\overline{n}\|}^{1},
Z\_{x:\overline{n}\|}^{1\text{(pure endow)}}) =
-A\_{x:\overline{n}\|}^{1} \cdot {}\_nE_x\\.

## Usage

``` r
cov_term_endow(x, n, i, tbl = NULL, model = NULL, ...)
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

Numeric vector of covariances.

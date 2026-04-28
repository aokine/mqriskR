# Covariance of term and deferred insurance PVs

Computes \\\mathrm{Cov}(Z\_{x:\overline{n}\|}^{1}, {}\_{n\mid}Z_x) =
-A\_{x:\overline{n}\|}^{1} \cdot {}\_{n\mid}A_x\\.

## Usage

``` r
cov_term_deferred(x, n, i, tbl = NULL, model = NULL, ...)
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

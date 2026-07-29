# Annuity-insurance relationships

Identities linking annual and continuous annuity functions to the
corresponding insurance functions.

Computes \\a_x = (v - A_x)/d\\.

Computes \\\ddot{a}\_x = (1 - A_x)/d\\.

Computes \\\bar{a}\_x = (1 - \bar{A}\_x)/\delta\\.

Computes \\a\_{x:\overline{n}\|} = (1 - A\_{x:\overline{n}\|})/d - 1 +
{}\_nE_x\\.

Computes \\\ddot{a}\_{x:\overline{n}\|} = (1 -
A\_{x:\overline{n}\|})/d\\.

Computes \\\bar{a}\_{x:\overline{n}\|} = (1 -
\bar{A}\_{x:\overline{n}\|})/\delta\\.

Computes \\{}\_{n\|}a_x = {}\_nE_x a\_{x+n}\\.

Computes \\{}\_{n\|}\ddot{a}\_x = {}\_nE_x \ddot{a}\_{x+n}\\.

Computes \\{}\_{n\|}\bar{a}\_x = {}\_nE_x \bar{a}\_{x+n}\\.

## Usage

``` r
annuity_identity_ax(x, i, model = NULL, ..., tbl = NULL)

annuity_identity_adotx(x, i, model = NULL, ..., tbl = NULL)

annuity_identity_abarx(x, i, model, ...)

annuity_identity_axn(x, n, i, model = NULL, ..., tbl = NULL)

annuity_identity_adotxn(x, n, i, model = NULL, ..., tbl = NULL)

annuity_identity_abarxn(x, n, i, model, ...)

annuity_identity_nax(
  x,
  n,
  i,
  model = NULL,
  ...,
  tbl = NULL,
  k_max = 5000,
  tol = 1e-12
)

annuity_identity_nadotx(
  x,
  n,
  i,
  model = NULL,
  ...,
  tbl = NULL,
  k_max = 5000,
  tol = 1e-12
)

annuity_identity_nabarx(x, n, i, model, ..., tol = 1e-10)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- model:

  Optional survival model name.

- ...:

  Additional model parameters.

- tbl:

  Optional life table object for discrete identities.

- n:

  Term or deferral period in years.

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

## Value

Numeric vector containing the annuity value computed from the
corresponding annuity-insurance identity.

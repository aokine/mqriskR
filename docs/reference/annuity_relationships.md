# Annuity-insurance relationships (Chapter 8)

This file provides the core Chapter 8 identities linking annual and
continuous annuity functions to the corresponding insurance functions.

Computes \\a_x = (v - A_x)/d\\.

Computes \\\ddot{a}\_x = (1 - A_x)/d\\.

Computes \\\bar{a}\_x = (1 - \bar{A}\_x)/\delta\\.

Computes \\a\_{x:\overline{n}\|} = \ddot{a}\_{x:\overline{n}\|} - 1 +
{}\_nE_x\\ together with \\\ddot{a}\_{x:\overline{n}\|} = (1 -
A\_{x:\overline{n}\|})/d\\.

Computes \\\ddot{a}\_{x:\overline{n}\|} = (1 -
A\_{x:\overline{n}\|})/d\\.

Computes \\\bar{a}\_{x:\overline{n}\|} = (1 -
\bar{A}\_{x:\overline{n}\|})/\delta\\.

Computes \\{}\_{n\|}a_x = {}\_nE_x\\ a\_{x+n}\\.

Computes \\{}\_{n\|}\ddot{a}\_x = {}\_nE_x\\ \ddot{a}\_{x+n}\\.

Computes \\{}\_{n\|}\bar{a}\_x = {}\_nE_x\\ \bar{a}\_{x+n}\\.

## Usage

``` r
annuity_identity_ax(x, i, model, ...)

annuity_identity_adotx(x, i, model, ...)

annuity_identity_abarx(x, i, model, ...)

annuity_identity_axn(x, n, i, model, ...)

annuity_identity_adotxn(x, n, i, model, ...)

annuity_identity_abarxn(x, n, i, model, ...)

annuity_identity_nax(x, n, i, model, ..., k_max = 5000, tol = 1e-12)

annuity_identity_nadotx(x, n, i, model, ..., k_max = 5000, tol = 1e-12)

annuity_identity_nabarx(x, n, i, model, ..., tol = 1e-10)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- model:

  Survival model.

- ...:

  Additional model parameters passed to the survival model.

- n:

  Term in years.

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

## Value

Numeric vector.

## Details

Included identities:

- whole life immediate: \\a_x = (v - A_x)/d\\

- whole life due: \\\ddot{a}\_x = (1 - A_x)/d\\

- whole life continuous: \\\bar{a}\_x = (1 - \bar{A}\_x)/\delta\\

- temporary immediate: \\a\_{x:\overline{n}\|} = (1 -
  A\_{x:\overline{n}\|})/d - 1 + {}\_nE_x\\

- temporary due: \\\ddot{a}\_{x:\overline{n}\|} = (1 -
  A\_{x:\overline{n}\|})/d\\

- temporary continuous: \\\bar{a}\_{x:\overline{n}\|} = (1 -
  \bar{A}\_{x:\overline{n}\|})/\delta\\

- deferred immediate: \\{}\_{n\mid}a_x = {}\_nE_x a\_{x+n}\\

- deferred due: \\{}\_{n\mid}\ddot{a}\_x = {}\_nE_x \ddot{a}\_{x+n}\\

- deferred continuous: \\{}\_{n\mid}\bar{a}\_x = {}\_nE_x
  \bar{a}\_{x+n}\\

These are wrapper functions that evaluate the Chapter 8 relationships
using the Chapter 7 insurance functions already implemented in the
package.

Hence \\a\_{x:\overline{n}\|} = (1 - A\_{x:\overline{n}\|})/d - 1 +
{}\_nE_x\\.

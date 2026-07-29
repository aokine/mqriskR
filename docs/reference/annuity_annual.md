# Annual annuity functions

Annual whole life, temporary, deferred, and actuarial accumulated value
annuity functions in immediate, due, and continuous forms.

Computes \\a_x = \sum\_{t=1}^{\infty} v^t {}\_t p_x\\.

Computes \\\ddot{a}\_x = \sum\_{t=0}^{\infty} v^t {}\_t p_x\\.

Computes \\\bar{a}\_x = \int_0^{\infty} v^t {}\_t p_x dt\\.

Computes \\a\_{x:\overline{n}\|} = \sum\_{t=1}^{n} v^t {}\_t p_x\\.

Computes \\\ddot{a}\_{x:\overline{n}\|} = \sum\_{t=0}^{n-1} v^t {}\_t
p_x\\.

Computes \\\bar{a}\_{x:\overline{n}\|} = \int_0^n v^t {}\_t p_x dt\\.

Computes \\{}\_{n\|}a_x = {}\_nE_x a\_{x+n}\\.

Computes \\{}\_{n\|}\ddot{a}\_x = {}\_nE_x \ddot{a}\_{x+n}\\.

Computes \\{}\_{n\|}\bar{a}\_x = {}\_nE_x \bar{a}\_{x+n}\\.

Computes \\s\_{x:\overline{n}\|} = a\_{x:\overline{n}\|} / {}\_nE_x\\.

Computes \\\ddot{s}\_{x:\overline{n}\|} = \ddot{a}\_{x:\overline{n}\|} /
{}\_nE_x\\.

Computes \\\bar{s}\_{x:\overline{n}\|} = \bar{a}\_{x:\overline{n}\|} /
{}\_nE_x\\.

## Usage

``` r
ax(x, i, model = NULL, ..., tbl = NULL, k_max = 5000, tol = 1e-12)

adotx(x, i, model = NULL, ..., tbl = NULL, k_max = 5000, tol = 1e-12)

abarx(x, i, model, ..., tol = 1e-10)

axn(x, n, i, model = NULL, ..., tbl = NULL)

adotxn(x, n, i, model = NULL, ..., tbl = NULL)

abarxn(x, n, i, model, ...)

nax(x, n, i, model = NULL, ..., tbl = NULL, k_max = 5000, tol = 1e-12)

nadotx(x, n, i, model = NULL, ..., tbl = NULL, k_max = 5000, tol = 1e-12)

nabarx(x, n, i, model, ..., tol = 1e-10)

sxn(x, n, i, model = NULL, ..., tbl = NULL)

sdotxn(x, n, i, model = NULL, ..., tbl = NULL)

sbarxn(x, n, i, model, ...)
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

  Optional life table object for annual discrete annuity functions.

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

- n:

  Term in years.

## Value

Numeric vector containing the requested annuity present value or
actuarial accumulated value.

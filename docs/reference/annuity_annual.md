# Annual annuity functions (Chapter 8)

Annual whole life, temporary, deferred, and actuarial accumulated value
annuity functions in immediate, due, and continuous forms.

Computes the annual whole life annuity-immediate \\a_x =
\sum\_{t=1}^{\infty} v^t \\ {}\_t p_x\\.

Computes the annual whole life annuity-due \\\ddot{a}\_x =
\sum\_{t=0}^{\infty} v^t \\ {}\_t p_x = 1 + a_x\\.

Computes the continuous whole life annuity \\\bar{a}\_x =
\int_0^{\infty} v^t \\ {}\_t p_x \\ dt\\.

Computes the annual temporary annuity-immediate \\a\_{x:\overline{n}\|}
= \sum\_{t=1}^{n} v^t \\ {}\_t p_x\\.

Computes the annual temporary annuity-due \\\ddot{a}\_{x:\overline{n}\|}
= \sum\_{t=0}^{n-1} v^t \\ {}\_t p_x\\.

Computes the continuous temporary annuity \\\bar{a}\_{x:\overline{n}\|}
= \int_0^n v^t \\ {}\_t p_x \\ dt\\.

Computes the annual deferred whole life annuity-immediate \\{}\_{n\|}a_x
= {}\_nE_x \\ a\_{x+n}\\.

Computes the annual deferred whole life annuity-due
\\{}\_{n\|}\ddot{a}\_x = {}\_nE_x \\ \ddot{a}\_{x+n}\\.

Computes the continuous deferred whole life annuity
\\{}\_{n\|}\bar{a}\_x = {}\_nE_x \\ \bar{a}\_{x+n}\\.

Computes \\s\_{x:\overline{n}\|} = a\_{x:\overline{n}\|} / {}\_nE_x\\.

Computes \\\ddot{s}\_{x:\overline{n}\|} = \ddot{a}\_{x:\overline{n}\|} /
{}\_nE_x\\.

Computes \\\bar{s}\_{x:\overline{n}\|} = \bar{a}\_{x:\overline{n}\|} /
{}\_nE_x\\.

## Usage

``` r
ax(x, i, model, ..., k_max = 5000, tol = 1e-12)

adotx(x, i, model, ..., k_max = 5000, tol = 1e-12)

abarx(x, i, model, ..., tol = 1e-10)

axn(x, n, i, model, ...)

adotxn(x, n, i, model, ...)

abarxn(x, n, i, model, ...)

nax(x, n, i, model, ..., k_max = 5000, tol = 1e-12)

nadotx(x, n, i, model, ..., k_max = 5000, tol = 1e-12)

nabarx(x, n, i, model, ..., tol = 1e-10)

sxn(x, n, i, model, ...)

sdotxn(x, n, i, model, ...)

sbarxn(x, n, i, model, ...)
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

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

- n:

  Term in years.

## Value

Numeric vector.

## Details

Naming convention follows Chapter 8 notation:

- `ax()` = \\a_x\\

- `adotx()` = \\\ddot{a}\_x\\

- `abarx()` = \\\bar{a}\_x\\

- `axn()` = \\a\_{x:\overline{n}\|}\\

- `adotxn()` = \\\ddot{a}\_{x:\overline{n}\|}\\

- `abarxn()` = \\\bar{a}\_{x:\overline{n}\|}\\

- `nax()` = \\{}\_{n\|}a_x\\

- `nadotx()` = \\{}\_{n\|}\ddot{a}\_x\\

- `nabarx()` = \\{}\_{n\|}\bar{a}\_x\\

- `sxn()` = \\s\_{x:\overline{n}\|}\\

- `sdotxn()` = \\\ddot{s}\_{x:\overline{n}\|}\\

- `sbarxn()` = \\\bar{s}\_{x:\overline{n}\|}\\

These functions work directly from the Chapter 5 survival model
functions and the Chapter 7 pure endowment function [`nEx()`](nEx.md).

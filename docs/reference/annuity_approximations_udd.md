# UDD annuity approximations

UDD-based approximations for Chapter 8 annuity functions.

Computes \$\$\ddot{a}\_x^{(m)} \approx \alpha(m)\ddot{a}\_x -
\beta(m).\$\$

Computes \$\$\ddot{a}\_{x:\overline{n}\|}^{(m)} \approx
\alpha(m)\ddot{a}\_{x:\overline{n}\|} - \beta(m)(1-{}\_nE_x).\$\$

Computes \$\${}\_{n\mid}\ddot{a}\_x^{(m)} \approx
\alpha(m)\\{}\_{n\mid}\ddot{a}\_x - \beta(m)\\{}\_nE_x.\$\$

Computes \$\$a_x^{(m)} \approx \alpha(m)a_x + \gamma(m).\$\$

Computes \$\$a\_{x:\overline{n}\|}^{(m)} \approx
\alpha(m)a\_{x:\overline{n}\|} + \gamma(m)(1-{}\_nE_x).\$\$

Computes \$\${}\_{n\mid}a_x^{(m)} \approx \alpha(m)\\{}\_{n\mid}a_x +
\gamma(m)\\{}\_nE_x.\$\$

Computes \$\$\ddot{s}\_{x:\overline{n}\|}^{(m)} \approx
\alpha(m)\ddot{s}\_{x:\overline{n}\|} -
\beta(m)\left(\frac{1}{{}\_nE_x}-1\right).\$\$

Computes \$\$s\_{x:\overline{n}\|}^{(m)} \approx
\alpha(m)s\_{x:\overline{n}\|} +
\gamma(m)\left(\frac{1}{{}\_nE_x}-1\right).\$\$

Computes \$\$\bar{a}\_x \approx \frac{id}{\delta^2}\ddot{a}\_x -
\frac{i-\delta}{\delta^2}.\$\$

Uses the identity \$\$\bar{a}\_{x:\overline{n}\|} \approx
\frac{1-\bar{A}\_{x:\overline{n}\|}}{\delta}\$\$ together with the
package's existing Chapter 7 UDD insurance approximation for
\\\bar{A}\_{x:\overline{n}\|}\\.

Computes \$\${}\_{n\mid}\bar{a}\_x \approx {}\_nE_x \\
\bar{a}\_{x+n}.\$\$

## Usage

``` r
adotx_m_udd(x, m, i, model, ...)

adotxn_m_udd(x, n, m, i, model, ...)

nadotx_m_udd(x, n, m, i, model, ...)

ax_m_udd(x, m, i, model, ...)

axn_m_udd(x, n, m, i, model, ...)

nax_m_udd(x, n, m, i, model, ...)

sdotxn_m_udd(x, n, m, i, model, ...)

sxn_m_udd(x, n, m, i, model, ...)

abarx_udd(x, i, model, ...)

abarxn_udd(x, n, i, model, ...)

nabarx_udd(x, n, i, model, ...)
```

## Arguments

- x:

  Age.

- m:

  Number of payments per year.

- i:

  Effective annual interest rate.

- model:

  Survival model.

- ...:

  Additional model parameters.

- n:

  Term.

## Value

Numeric vector.

## Details

These functions implement the standard Uniform Distribution of Deaths
approximations linking annual, m-thly, and continuous annuity values.

The exported functions documented on this page are:

- `ax_m_udd()`

- `axn_m_udd()`

- `nax_m_udd()`

- `adotx_m_udd()`

- `adotxn_m_udd()`

- `nadotx_m_udd()`

- `sxn_m_udd()`

- `sdotxn_m_udd()`

- `abarx_udd()`

- `abarxn_udd()`

- `nabarx_udd()`

Note that this function relies on the already-existing
[`Abarxn_udd()`](Abarxn_udd.md) implementation in the package, so extra
survival-model arguments are not used.

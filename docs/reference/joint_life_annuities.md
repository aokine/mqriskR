# Joint-life annuities

Computes temporary and whole-life joint-life annuities-due and
annuities-immediate for two independent lives.

## Usage

``` r
adotxyn(x, y, n, i, tbl = NULL, model = NULL, ...)

axyn(x, y, n, i, tbl = NULL, model = NULL, ...)

adotxy(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)

axy(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000L, tol = 1e-12)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

- n:

  Nonnegative integer term. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the survival model.

- k_max:

  Maximum summation horizon used for non-terminating parametric survival
  models.

- tol:

  Positive convergence tolerance for whole-life summations.

## Value

A numeric vector of annuity actuarial present values.

## Details

`adotxyn()` computes an \\n\\-year temporary joint-life annuity-due.

`axyn()` computes an \\n\\-year temporary joint-life annuity-immediate.

`adotxy()` computes a whole-life joint-life annuity-due.

`axy()` computes a whole-life joint-life annuity-immediate.

All calculations assume independence between the two future lifetimes.

The temporary annuity-due is \$\$ \ddot{a}\_{xy:\overline{n}\|} =
\sum\_{k=0}^{n-1} v^k\\{}\_kp\_{xy}. \$\$

The temporary annuity-immediate is \$\$ a\_{xy:\overline{n}\|} =
\sum\_{k=1}^{n} v^k\\{}\_kp\_{xy}. \$\$

For life-table calculations, the whole-life sums terminate at the last
duration supported jointly by the two lives. For parametric models, the
sums are evaluated until convergence or until `k_max` is reached.

## Examples

``` r
adotxyn(
  x = 40,
  y = 50,
  n = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 6.956659

axyn(
  x = 40,
  y = 50,
  n = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 6.365935

adotxy(
  x = 40,
  y = 50,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 10.96154

axy(
  x = 40,
  y = 50,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 9.961536
```

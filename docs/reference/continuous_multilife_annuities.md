# Continuous multi-life annuities

Computes continuous joint-life, last-survivor, and reversionary
whole-life annuities for two independent lives.

## Usage

``` r
abarxy(x, y, i, tbl = NULL, model = NULL, ...)

abarxybar(x, y, i, tbl = NULL, model = NULL, ...)

abarx_y(x, y, i, tbl = NULL, model = NULL, ...)

abary_x(x, y, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object retained for backward compatibility.
  Continuous calculations currently require `model`.

- model:

  Parametric survival model.

- ...:

  Additional parameters passed to the survival functions.

## Value

A numeric vector of actuarial present values.

## Details

`abarxy()` computes the joint-life annuity.

`abarxybar()` computes the last-survivor annuity.

`abarx_y()` computes the reversionary annuity payable to the second life
after the death of the first.

`abary_x()` computes the reversionary annuity payable to the first life
after the death of the second.

These functions require a parametric survival model.

Under independence, `abarxy()` represents the continuous joint-life
annuity, payable while both lives survive.

The last-survivor annuity satisfies \$\$\bar{a}\_{\overline{xy}} =
\bar{a}\_x + \bar{a}\_y - \bar{a}\_{xy}.\$\$

The reversionary annuities are obtained as the difference between the
corresponding single-life annuity and the joint-life annuity.

## Examples

``` r
abarxy(
  x = 40,
  y = 50,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 10.45444

abarxybar(
  x = 40,
  y = 50,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 16.24185
```

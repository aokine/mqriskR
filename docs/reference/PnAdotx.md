# Net premium for a deferred annuity-due

Computes \$\$ P({}\_{n\|}\ddot{a}\_x) = \frac{{}\_{n\|}\ddot{a}\_x}
{\ddot{a}\_{x:\overline{n}\|}}. \$\$

## Usage

``` r
PnAdotx(x, n, i, model = NULL, ..., tbl = NULL)
```

## Arguments

- x:

  Issue age. May be scalar or vector.

- n:

  Positive integer deferral period. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the actuarial functions.

- tbl:

  Optional life table object. Supply by name.

## Value

Numeric vector of net annual premiums.

## Examples

``` r
PnAdotx(
  40,
  n = 20,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.2651853
```

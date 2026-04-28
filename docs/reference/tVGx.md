# Whole life gross premium reserve

Computes the Chapter 11 prospective gross premium reserve for a fully
discrete whole life insurance with annual premiums and renewal expenses.

Computes the Chapter 11 prospective gross premium reserve for a fully
discrete whole life insurance with annual premiums and renewal expenses.

## Usage

``` r
tVGx(
  x,
  t,
  i,
  G,
  benefit = 1,
  renewal_premium_pct = 0,
  renewal_policy_exp = 0,
  settlement_exp = 0,
  tbl = NULL,
  model = NULL,
  ...
)

tVGx(
  x,
  t,
  i,
  G,
  benefit = 1,
  renewal_premium_pct = 0,
  renewal_policy_exp = 0,
  settlement_exp = 0,
  tbl = NULL,
  model = NULL,
  ...
)
```

## Arguments

- x:

  Issue age.

- t:

  Duration.

- i:

  Effective annual interest rate.

- G:

  Gross annual premium.

- benefit:

  Insurance amount. Default 1.

- renewal_premium_pct:

  Renewal percent-of-premium expense.

- renewal_policy_exp:

  Renewal per-policy expense.

- settlement_exp:

  Settlement expense at death.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

Numeric vector.

## Details

This function is intended for durations after issue, where future
expenses are modeled through renewal premium expenses, renewal
per-policy expenses, and settlement expense.

## Examples

``` r
tVGx(
  x = 40, t = 10, i = 0.05, G = 0.03,
  benefit = 1, renewal_premium_pct = 0.10,
  renewal_policy_exp = 0.002, settlement_exp = 0.02,
  model = "uniform", omega = 100
)
#> [1] 0.0391081
tVGx(
  x = 40, t = 10, i = 0.05, G = 0.03,
  benefit = 1, renewal_premium_pct = 0.10,
  renewal_policy_exp = 0.002, settlement_exp = 0.02,
  model = "uniform", omega = 100
)
#> [1] 0.0391081
```

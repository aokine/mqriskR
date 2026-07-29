# Whole life gross premium and expense reserves

Computes prospective gross premium and expense reserves for fully
discrete whole life insurance after issue.

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

tVEx(
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

  Issue age. May be scalar or vector.

- t:

  Nonnegative integer duration. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- G:

  Gross annual premium. May be scalar or vector.

- benefit:

  Insurance benefit amount.

- renewal_premium_pct:

  Renewal percent-of-premium expense in \\\[0,1\]\\.

- renewal_policy_exp:

  Renewal per-policy expense.

- settlement_exp:

  Settlement expense paid at death.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the actuarial functions.

## Value

A numeric vector of reserve values.

## Details

`tVGx()` computes the prospective gross premium reserve.

`tVEx()` computes the corresponding expense reserve, defined as the
difference between the gross premium reserve and the net benefit
reserve.

The gross premium reserve is calculated as

\$\$ {}\_tV_x^G = (b+s)A\_{x+t} -
\left\[(1-r)G-e\right\]\ddot{a}\_{x+t}, \$\$

where

- \\b\\ is the insurance benefit,

- \\s\\ is the settlement expense,

- \\G\\ is the gross annual premium,

- \\r\\ is the renewal percent-of-premium expense, and

- \\e\\ is the renewal per-policy expense.

The expense reserve is obtained as the gross premium reserve minus the
corresponding net benefit reserve.

## Examples

``` r
tVGx(
  x = 40,
  t = 10,
  i = 0.05,
  G = 0.03,
  benefit = 1,
  renewal_premium_pct = 0.10,
  renewal_policy_exp = 0.002,
  settlement_exp = 0.02,
  model = "uniform",
  omega = 100
)
#> [1] 0.0391081

tVEx(
  x = 40,
  t = 10,
  i = 0.05,
  G = 0.03,
  benefit = 1,
  renewal_premium_pct = 0.10,
  renewal_policy_exp = 0.002,
  settlement_exp = 0.02,
  model = "uniform",
  omega = 100
)
#> [1] -0.03339664
```

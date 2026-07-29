# Total gross gain for a discrete insurance contract

Computes total gain during one policy year under gross premiums, gross
reserves, actual mortality, actual interest, and actual expenses.

## Usage

``` r
GTg_disc(
  VtG,
  Vt1G,
  G,
  i_actual,
  q_actual,
  r_actual = 0,
  e_actual = 0,
  s_actual = 0,
  b = 1
)
```

## Arguments

- VtG:

  Gross reserve at duration `t`.

- Vt1G:

  Gross reserve at duration `t + 1`.

- G:

  Gross premium.

- i_actual:

  Actual annual effective interest rate.

- q_actual:

  Actual mortality probability.

- r_actual:

  Actual percent-of-premium expense rate.

- e_actual:

  Actual per-policy expense.

- s_actual:

  Actual settlement expense.

- b:

  Benefit amount.

## Value

A numeric vector of total gains.

## Examples

``` r
GTg_disc(
  VtG = 0.10,
  Vt1G = 0.12,
  G = 0.02,
  i_actual = 0.05,
  q_actual = 0.01,
  r_actual = 0.03,
  e_actual = 0,
  s_actual = 0.01,
  b = 1
)
#> [1] -0.00353
```

# Total gross gain for a discrete insurance contract

Computes the Chapter 11 total gain under gross premiums and gross
reserves.

Computes the Chapter 11 total gain under gross premiums and gross
reserves.

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

  Gross reserve at duration t.

- Vt1G:

  Gross reserve at duration t+1.

- G:

  Gross premium.

- i_actual:

  Actual annual effective interest rate.

- q_actual:

  Actual mortality rate.

- r_actual:

  Actual percent-of-premium expense rate.

- e_actual:

  Actual per-policy expense.

- s_actual:

  Actual settlement expense.

- b:

  Benefit amount. Default 1.

## Value

Numeric vector.

Numeric vector.

## Examples

``` r
GTg_disc(
  VtG = 0.10, Vt1G = 0.12, G = 0.02,
  i_actual = 0.05, q_actual = 0.01,
  r_actual = 0.03, e_actual = 0, s_actual = 0.01, b = 1
)
#> [1] -0.00353
GTg_disc(
  VtG = 0.10, Vt1G = 0.12, G = 0.02,
  i_actual = 0.05, q_actual = 0.01,
  r_actual = 0.03, e_actual = 0, s_actual = 0.01, b = 1
)
#> [1] -0.00353
```

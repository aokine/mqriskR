# Ordered decomposition of gross gain

Decomposes gross gain into interest, mortality, and expense components
in a user-specified order.

## Usage

``` r
decompGg_disc(
  VtG,
  Vt1G,
  G,
  i_assumed,
  q_assumed,
  r_assumed = 0,
  e_assumed = 0,
  s_assumed = 0,
  i_actual,
  q_actual,
  r_actual = 0,
  e_actual = 0,
  s_actual = 0,
  b = 1,
  order = c("interest", "mortality", "expense")
)
```

## Arguments

- VtG:

  Gross reserve at duration `t`.

- Vt1G:

  Gross reserve at duration `t + 1`.

- G:

  Gross premium.

- i_assumed:

  Assumed annual effective interest rate.

- q_assumed:

  Assumed mortality probability.

- r_assumed:

  Assumed percent-of-premium expense rate.

- e_assumed:

  Assumed per-policy expense.

- s_assumed:

  Assumed settlement expense.

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

- order:

  Character vector containing `"interest"`, `"mortality"`, and
  `"expense"` exactly once.

## Value

For scalar input, a named numeric vector. For vectorized input, a
numeric matrix with columns `total_gain`, `interest`, `mortality`,
`expense`, and `check`.

## Details

For scalar input, the function returns a named numeric vector. For
vectorized input, it returns a numeric matrix with one row per
calculation.

## Examples

``` r
decompGg_disc(
  VtG = 3950.73,
  Vt1G = 4607.07,
  G = 685,
  i_assumed = 0.06,
  q_assumed = 0.00592,
  r_assumed = 0.05,
  e_assumed = 0,
  s_assumed = 300,
  i_actual = 0.065,
  q_actual = 0.005,
  r_actual = 0.06,
  e_actual = 0,
  s_actual = 100,
  b = 50000
)
#> total_gain   interest  mortality    expense      check 
#>   58.74630   23.00405   42.03750   -6.29525   58.74630 
```

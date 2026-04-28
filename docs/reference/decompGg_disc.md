# Ordered gross gain decomposition

Decomposes total gross gain into interest, mortality, and expense
components in a user-specified order.

Decomposes total gross gain into interest, mortality, and expense
components in a user-specified order. The first two components are
computed sequentially, and the last component is taken as the balancing
item so that the components sum exactly to total gain.

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

  Gross reserve at duration t.

- Vt1G:

  Gross reserve at duration t+1.

- G:

  Gross premium.

- i_assumed:

  Assumed annual effective interest rate.

- q_assumed:

  Assumed mortality rate.

- r_assumed:

  Assumed percent-of-premium expense rate.

- e_assumed:

  Assumed per-policy expense.

- s_assumed:

  Assumed settlement expense.

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

- order:

  Character vector giving the order of decomposition.

## Value

Named numeric vector.

Named numeric vector.

## Examples

``` r
decompGg_disc(
  VtG = 3950.73, Vt1G = 4607.07, G = 685,
  i_assumed = 0.06, q_assumed = 0.00592,
  r_assumed = 0.05, e_assumed = 0, s_assumed = 300,
  i_actual = 0.065, q_actual = 0.005,
  r_actual = 0.06, e_actual = 0, s_actual = 100,
  b = 50000,
  order = c("interest", "mortality", "expense")
)
#> total_gain   interest  mortality    expense      check 
#>   58.74630   23.00405   42.03750   -6.29525   58.74630 
decompGg_disc(
  VtG = 3950.73, Vt1G = 4607.07, G = 685,
  i_assumed = 0.06, q_assumed = 0.00592,
  r_assumed = 0.05, e_assumed = 0, s_assumed = 300,
  i_actual = 0.065, q_actual = 0.005,
  r_actual = 0.06, e_actual = 0, s_actual = 100,
  b = 50000,
  order = c("interest", "mortality", "expense")
)
#> total_gain   interest  mortality    expense      check 
#>   58.74630   23.00405   42.03750   -6.29525   58.74630 
```

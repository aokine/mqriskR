# Present value of deterministic cash flows using spot rates

Discounts each cash flow using the spot rate matched to its payment
time. Time-zero cash flows are not discounted.

## Usage

``` r
pv_spot_cashflows(
  amounts,
  times,
  spot,
  compounding = c("annual", "semiannual")
)
```

## Arguments

- amounts:

  Numeric vector of cash-flow amounts.

- times:

  Numeric vector of payment times in years.

- spot:

  Numeric vector of spot rates matched elementwise to `times`. The value
  corresponding to a time-zero cash flow is ignored.

- compounding:

  Character string equal to `"annual"` or `"semiannual"`.

## Value

A numeric scalar.

## Examples

``` r
pv_spot_cashflows(
  amounts = c(200000, 50000, 50000, 100000),
  times = c(0, 0.5, 1, 2),
  spot = c(0, 0.02440, 0.02601, 0.02936),
  compounding = "semiannual"
)
#> [1] 392459.1
```

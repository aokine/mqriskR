# Present value of cash flows using spot rates

Discounts deterministic cash flows using spot rates matched to their
maturities.

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

  Numeric vector of cash flow amounts.

- times:

  Numeric vector of payment times in years.

- spot:

  Numeric vector of spot rates matched elementwise to `times`. Use 0 for
  any time-0 entry.

- compounding:

  Either `"annual"` or `"semiannual"`.

## Value

A numeric scalar.

## Details

For annual compounding, each positive-time cash flow at time \\t\\ is
discounted by \\(1+z_t)^{-t}\\.

For semiannual nominal compounding, each positive-time cash flow at time
\\t\\ is discounted by \\(1+z_t/2)^{-2t}\\.

Time-0 cash flows are left undiscounted.

## Examples

``` r
pv_spot_cashflows(
  amounts = c(200000, 50000, 50000, 100000),
  times   = c(0, 0.5, 1, 2),
  spot    = c(0, 0.02440, 0.02601, 0.02936),
  compounding = "semiannual"
)
#> [1] 392459.1
```

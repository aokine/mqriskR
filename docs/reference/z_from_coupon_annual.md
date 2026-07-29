# Bootstrap annual effective spot rates

Bootstraps annual effective spot rates from par coupon yields at
consecutive integer maturities.

## Usage

``` r
z_from_coupon_annual(maturity, coupon_yield, par = 1000)
```

## Arguments

- maturity:

  Numeric vector of positive integer maturities in strictly increasing
  order. Maturities must be consecutive and begin at 1.

- coupon_yield:

  Numeric vector of annual effective par coupon yields. Values must be
  greater than `-1`.

- par:

  Positive scalar par value.

## Value

A numeric vector of annual effective spot rates.

## Examples

``` r
maturity <- 1:4
coupon_yield <- c(0.02, 0.04, 0.06, 0.08)
z_from_coupon_annual(maturity, coupon_yield)
#> [1] 0.02000000 0.04040808 0.06169260 0.08447397
```

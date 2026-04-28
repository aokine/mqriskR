# Bootstrap annual spot rates from annual coupon-bond yields

Bootstraps annual effective zero-coupon yields from par annual
coupon-bearing bond yields of the same maturities.

## Usage

``` r
z_from_coupon_annual(maturity, coupon_yield, par = 1000)
```

## Arguments

- maturity:

  Integer vector of maturities in years, in increasing order.

- coupon_yield:

  Numeric vector of annual coupon yields.

- par:

  Par value of each bond.

## Value

Numeric vector of annual effective spot rates.

## Examples

``` r
maturity <- 1:4
coupon_yield <- c(0.02, 0.04, 0.06, 0.08)
z_from_coupon_annual(maturity, coupon_yield)
#> [1] 0.02000000 0.04040808 0.06169260 0.08447397
```

# Bootstrap semiannual nominal spot rates from coupon-bond yields

Bootstraps the semiannual nominal annual zero-coupon yields from par
coupon-bearing bond yields of the same maturities.

## Usage

``` r
z_from_coupon_semi(maturity, coupon_yield, par = 1000)
```

## Arguments

- maturity:

  Numeric vector of maturities in years, typically `0.5, 1.0, 1.5, ...`,
  in increasing order.

- coupon_yield:

  Numeric vector of nominal annual coupon yields convertible
  semiannually.

- par:

  Par value of each bond.

## Value

Numeric vector of semiannual nominal annual spot rates.

## Details

Both coupon yields and spot yields are interpreted as nominal annual
rates convertible semiannually.

## Examples

``` r
maturity <- c(0.5, 1.0, 1.5, 2.0)
coupon_yield <- c(0.0244, 0.0260, 0.0276, 0.0293)
z_from_coupon_semi(maturity, coupon_yield)
#> [1] 0.02440000 0.02601041 0.02762959 0.02936142
```

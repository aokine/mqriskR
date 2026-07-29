# Bootstrap semiannual nominal spot rates

Bootstraps nominal annual spot rates convertible semiannually from par
coupon yields at consecutive half-year maturities.

## Usage

``` r
z_from_coupon_semi(maturity, coupon_yield, par = 1000)
```

## Arguments

- maturity:

  Numeric vector of positive maturities in years, in strictly increasing
  order. Maturities must be consecutive multiples of `0.5`.

- coupon_yield:

  Numeric vector of nominal annual par coupon yields convertible
  semiannually. Values must be greater than `-2`.

- par:

  Positive scalar par value.

## Value

A numeric vector of nominal annual spot rates convertible semiannually.

## Examples

``` r
maturity <- c(0.5, 1.0, 1.5, 2.0)
coupon_yield <- c(0.0244, 0.0260, 0.0276, 0.0293)
z_from_coupon_semi(maturity, coupon_yield)
#> [1] 0.02440000 0.02601041 0.02762959 0.02936142
```

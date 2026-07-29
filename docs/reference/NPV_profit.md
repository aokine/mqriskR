# Net present value of a profit signature

Discounts a profit signature at one or more annual effective risk
discount rates.

## Usage

``` r
NPV_profit(Pi, r)
```

## Arguments

- Pi:

  Numeric profit-signature vector.

- r:

  Annual effective risk discount rate. May be scalar or vector; values
  must be greater than `-1`.

## Value

A numeric vector with one value for each rate in `r`.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
NPV_profit(Pi, r = 0.10)
#> [1] 6.034711
NPV_profit(Pi, r = c(0.08, 0.10, 0.12))
#> [1] 6.800450 6.034711 5.313388
```

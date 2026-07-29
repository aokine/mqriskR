# Discounted payback period

Returns the first duration at which cumulative discounted profit is
nonnegative.

## Usage

``` r
discounted_payback_period(Pi, r)
```

## Arguments

- Pi:

  Numeric profit-signature vector.

- r:

  Annual effective risk discount rate. May be scalar or vector; values
  must be greater than `-1`.

## Value

An integer vector with one value for each discount rate. An element is
`NA_integer_` when payback is not reached.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
discounted_payback_period(Pi, r = 0.10)
#> [1] 3
```

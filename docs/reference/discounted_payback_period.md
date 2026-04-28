# Discounted payback period

Returns the first duration \\t\\ for which the partial net present value
\\NPV(t)\\ is nonnegative.

## Usage

``` r
discounted_payback_period(Pi, r)
```

## Arguments

- Pi:

  Profit signature vector.

- r:

  Risk discount rate.

## Value

Integer scalar, or `NA_integer_` if the payback period is not reached.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
discounted_payback_period(Pi, r = 0.10)
#> [1] 3
```

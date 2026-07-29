# Internal rate of return

Computes a root of the profit-signature net present value.

## Usage

``` r
IRR_profit(Pi, interval = c(0, 1), tol = .Machine$double.eps^0.5)
```

## Arguments

- Pi:

  Numeric profit-signature vector.

- interval:

  Numeric vector of length two giving the root-search interval. Its
  lower endpoint must be greater than `-1`.

- tol:

  Positive scalar tolerance passed to
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html).

## Value

A numeric scalar.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
IRR_profit(Pi)
#> [1] 0.3163508
```

# Internal rate of return of a profit signature

Computes the internal rate of return (IRR), defined as the rate \\r\\
for which the net present value is zero.

## Usage

``` r
IRR_profit(Pi, interval = c(0, 1), tol = .Machine$double.eps^0.5)
```

## Arguments

- Pi:

  Profit signature vector.

- interval:

  Numeric vector of length 2 giving the search interval for
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html).

- tol:

  Tolerance passed to
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html).

## Value

Numeric scalar.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
IRR_profit(Pi)
#> [1] 0.3163508
```

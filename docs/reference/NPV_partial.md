# Partial net present values

Computes cumulative discounted profits through each duration. A scalar
discount rate returns a named vector; vectorized rates return a matrix.

## Usage

``` r
NPV_partial(Pi, r)
```

## Arguments

- Pi:

  Numeric profit-signature vector.

- r:

  Annual effective risk discount rate. May be scalar or vector; values
  must be greater than `-1`.

## Value

A named numeric vector for scalar `r`, or a numeric matrix for
vectorized `r`.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
NPV_partial(Pi, r = 0.10)
#>      NPV(0)      NPV(1)      NPV(2)      NPV(3) 
#> -15.0000000  -7.3454545  -0.4115702   6.0347107 
```

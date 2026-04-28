# Partial net present values

Computes the sequence \\NPV(0), NPV(1), \dots, NPV(n)\\ of partial net
present values from a profit signature.

## Usage

``` r
NPV_partial(Pi, r)
```

## Arguments

- Pi:

  Profit signature vector.

- r:

  Risk discount rate.

## Value

Numeric vector.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
NPV_partial(Pi, r = 0.10)
#>      NPV(0)      NPV(1)      NPV(2)      NPV(3) 
#> -15.0000000  -7.3454545  -0.4115702   6.0347107 
```

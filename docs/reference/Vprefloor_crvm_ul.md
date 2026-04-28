# Pre-floor CRVM reserve for universal life

Computes the pre-floor CRVM reserve from Equation (16.17).

## Usage

``` r
Vprefloor_crvm_ul(r, pvfb_minus_pvfp)
```

## Arguments

- r:

  Ratio \\r_t\\.

- pvfb_minus_pvfp:

  Difference \\(PVFB)\_t - (PVFP)\_t\\.

## Value

Numeric scalar.

## Examples

``` r
Vprefloor_crvm_ul(r = 0.33506, pvfb_minus_pvfp = 70)
#> [1] 23.4542
```

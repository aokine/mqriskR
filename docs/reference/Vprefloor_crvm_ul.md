# Pre-floor CRVM reserve

Computes the pre-floor reserve by multiplying the funding ratio by the
difference between the present value of future benefits and future
premiums.

## Usage

``` r
Vprefloor_crvm_ul(r, pvfb_minus_pvfp)
```

## Arguments

- r:

  Funding ratio. Values must lie in `[0, 1]`.

- pvfb_minus_pvfp:

  Numeric difference between the present value of future benefits and
  future premiums.

## Value

A numeric vector.

## Examples

``` r
Vprefloor_crvm_ul(r = 0.33506, pvfb_minus_pvfp = 70)
#> [1] 23.4542
```

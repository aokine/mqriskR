# Guaranteed maturity fund roll-forward

Computes a one-period guaranteed maturity fund roll-forward. Arguments
may be scalars or vectors and follow common-length recycling.

## Usage

``` r
GMF_rollforward_ul(GMF_prev, GMP, r, policy_charge, i)
```

## Arguments

- GMF_prev:

  Prior guaranteed maturity fund.

- GMP:

  Guaranteed maturity premium.

- r:

  Percent-of-premium expense rate in `[0, 1]`.

- policy_charge:

  Guaranteed policy charge.

- i:

  Guaranteed annual effective interest rate greater than `-1`.

## Value

A numeric vector.

## Examples

``` r
GMF_rollforward_ul(140.40, 14.49, 0.04, 11.80, 0.03)
#> [1] 146.7857
```

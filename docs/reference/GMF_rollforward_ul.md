# Guaranteed maturity fund roll-forward

Computes the one-period guaranteed maturity fund roll-forward used in
Example 16.9.

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

  Expense factor applied to GMP.

- policy_charge:

  Guaranteed policy charge.

- i:

  Guaranteed interest rate.

## Value

Numeric scalar.

## Examples

``` r
GMF_rollforward_ul(140.40, 14.49, 0.04, 11.80, 0.03)
#> [1] 146.7857
```

# Entry Age Normal normal cost

Computes Entry Age Normal normal cost as total benefit APV divided by an
active-service annuity-due factor.

## Usage

``` r
NC_EAN_db(APV_total, adue_active)
```

## Arguments

- APV_total:

  Nonnegative total actuarial present value of benefits.

- adue_active:

  Positive active-service annuity-due factor.

## Value

A numeric vector.

## Examples

``` r
NC_EAN_db(APV_total = 25000, adue_active = 15)
#> [1] 1666.667
```

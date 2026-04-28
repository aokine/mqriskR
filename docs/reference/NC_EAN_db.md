# Entry Age Normal normal cost for a DB plan

Computes Equation (18.13).

## Usage

``` r
NC_EAN_db(APV_total, adue_active)
```

## Arguments

- APV_total:

  Total actuarial present value of benefits.

- adue_active:

  Active-service annuity-due factor.

## Value

Entry Age Normal normal cost.

## Examples

``` r
NC_EAN_db(APV_total = 25000, adue_active = 15)
#> [1] 1666.667
```

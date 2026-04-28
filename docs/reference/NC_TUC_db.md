# Traditional Unit Credit normal cost for a DB plan

Computes the TUC normal cost as the APV of the current year's accrual.

## Usage

``` r
NC_TUC_db(accrual_benefit, v_to_ret, p_surv, adue_ret)
```

## Arguments

- accrual_benefit:

  Benefit accrued in the current year.

- v_to_ret:

  Discount factor to retirement.

- p_surv:

  Active-service survival probability to retirement.

- adue_ret:

  Retirement annuity factor.

## Value

TUC normal cost.

## Examples

``` r
NC_TUC_db(accrual_benefit = 1560, v_to_ret = 0.5, p_surv = 0.9, adue_ret = 12)
#> [1] 8424
```

# Traditional Unit Credit accrued liability for a DB plan

Computes the TUC accrued liability as the APV of the accrued benefit.

## Usage

``` r
AAL_TUC_db(accrued_benefit, v_to_ret, p_surv, adue_ret)
```

## Arguments

- accrued_benefit:

  Accrued benefit at the valuation date.

- v_to_ret:

  Discount factor to retirement.

- p_surv:

  Active-service survival probability to retirement.

- adue_ret:

  Retirement annuity factor.

## Value

TUC accrued liability.

## Examples

``` r
AAL_TUC_db(accrued_benefit = 12000, v_to_ret = 0.5, p_surv = 0.9, adue_ret = 12)
#> [1] 64800
```

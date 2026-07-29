# Traditional Unit Credit normal cost

Computes the actuarial present value of the benefit accrued during the
current year.

## Usage

``` r
NC_TUC_db(accrual_benefit, v_to_ret, p_surv, adue_ret)
```

## Arguments

- accrual_benefit:

  Nonnegative benefit accrued during the current year.

- v_to_ret:

  Nonnegative discount factor from the valuation date to retirement.

- p_surv:

  Survival or active-service probability to retirement in `[0, 1]`.

- adue_ret:

  Positive retirement annuity-due factor.

## Value

A numeric vector.

## Examples

``` r
NC_TUC_db(
  accrual_benefit = 1560,
  v_to_ret = 0.5,
  p_surv = 0.9,
  adue_ret = 12
)
#> [1] 8424
```

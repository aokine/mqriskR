# Traditional Unit Credit accrued liability

Computes the actuarial present value of the benefit accrued through the
valuation date.

## Usage

``` r
AAL_TUC_db(accrued_benefit, v_to_ret, p_surv, adue_ret)
```

## Arguments

- accrued_benefit:

  Nonnegative accrued benefit.

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
AAL_TUC_db(
  accrued_benefit = 12000,
  v_to_ret = 0.5,
  p_surv = 0.9,
  adue_ret = 12
)
#> [1] 64800
```

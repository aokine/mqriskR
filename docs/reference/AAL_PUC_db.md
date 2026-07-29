# Projected Unit Credit accrued liability

Computes the actuarial present value of the portion of the projected
benefit attributed to past service.

## Usage

``` r
AAL_PUC_db(
  projected_benefit,
  past_service,
  total_service,
  v_to_ret,
  p_surv,
  adue_ret
)
```

## Arguments

- projected_benefit:

  Nonnegative projected annual benefit at retirement.

- past_service:

  Nonnegative service completed through the valuation date.

- total_service:

  Positive total service at retirement.

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
AAL_PUC_db(
  projected_benefit = 30000,
  past_service = 10,
  total_service = 30,
  v_to_ret = 0.5,
  p_surv = 0.9,
  adue_ret = 12
)
#> [1] 54000
```

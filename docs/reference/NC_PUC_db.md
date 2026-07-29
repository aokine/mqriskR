# Projected Unit Credit normal cost

Computes the actuarial present value of the portion of the projected
benefit attributed to the current year of service.

## Usage

``` r
NC_PUC_db(projected_benefit, total_service, v_to_ret, p_surv, adue_ret)
```

## Arguments

- projected_benefit:

  Nonnegative projected annual benefit at retirement.

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
NC_PUC_db(
  projected_benefit = 30000,
  total_service = 30,
  v_to_ret = 0.5,
  p_surv = 0.9,
  adue_ret = 12
)
#> [1] 5400
```

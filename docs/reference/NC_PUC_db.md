# Projected Unit Credit normal cost for a DB plan

Computes the PUC normal cost as the APV of the portion of the projected
benefit attributed to the current year of service.

## Usage

``` r
NC_PUC_db(projected_benefit, total_service, v_to_ret, p_surv, adue_ret)
```

## Arguments

- projected_benefit:

  Projected benefit at retirement.

- total_service:

  Total service at retirement.

- v_to_ret:

  Discount factor to retirement.

- p_surv:

  Active-service survival probability to retirement.

- adue_ret:

  Retirement annuity factor.

## Value

PUC normal cost.

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

# Projected Unit Credit accrued liability for a DB plan

Computes the PUC accrued liability as the APV of the portion of the
projected benefit attributed to past service.

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

  Projected benefit at retirement.

- past_service:

  Past service completed.

- total_service:

  Total service at retirement.

- v_to_ret:

  Discount factor to retirement.

- p_surv:

  Active-service survival probability to retirement.

- adue_ret:

  Retirement annuity factor.

## Value

PUC accrued liability.

## Examples

``` r
AAL_PUC_db(projected_benefit = 30000, past_service = 10, total_service = 30,
v_to_ret = 0.5, p_surv = 0.9, adue_ret = 12)
#> [1] 54000
```

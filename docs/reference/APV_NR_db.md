# Actuarial present value of a normal retirement benefit

Computes the actuarial present value of a projected annual retirement
benefit. Scalar arguments are recycled to a common length.

## Usage

``` r
APV_NR_db(PABz, v_to_ret, p_surv, adue_ret)
```

## Arguments

- PABz:

  Nonnegative projected annual benefit at retirement.

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
APV_NR_db(
  PABz = 108008.66,
  v_to_ret = 1 / 1.06^30,
  p_surv = 0.8,
  adue_ret = 12
)
#> [1] 180531.9
```

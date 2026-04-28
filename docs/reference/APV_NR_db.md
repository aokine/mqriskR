# APV of normal retirement benefit for a DB plan

Computes Equation (18.10).

## Usage

``` r
APV_NR_db(PABz, v_to_ret, p_surv, adue_ret)
```

## Arguments

- PABz:

  Projected annual benefit at retirement.

- v_to_ret:

  Discount factor from current age to retirement.

- p_surv:

  Active-service survival probability to retirement.

- adue_ret:

  Retirement annuity factor.

## Value

Actuarial present value of the normal retirement benefit.

## Examples

``` r
APV_NR_db(PABz = 108008.66, v_to_ret = 1 / 1.06^30, p_surv = 0.8, adue_ret = 12)
#> [1] 180531.9
```

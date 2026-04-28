# Total gain for a continuous-style one-step recursion

Total gain for a continuous-style one-step recursion

## Usage

``` r
GT_cont(Vt, Vt1, P, delta_actual, p_actual, benefit = 0, h = 1)
```

## Arguments

- Vt:

  Reserve at time t.

- Vt1:

  Reserve at time t+h.

- P:

  Premium rate.

- delta_actual:

  Actual force of interest.

- p_actual:

  Actual survival probability over the step.

- benefit:

  Benefit paid at start of step. Default 0.

- h:

  Step length. Default 1.

## Value

Numeric vector.

## Examples

``` r
GT_cont(Vt = 10, Vt1 = 11, P = 1, delta_actual = 0.05, p_actual = 0.99)
#> [1] 0.6739821
```

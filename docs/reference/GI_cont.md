# Interest gain helper for continuous-style recursion

Interest gain helper for continuous-style recursion

## Usage

``` r
GI_cont(Vt, Vt1, P, delta_actual, p_assumed, benefit = 0, h = 1)
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

- p_assumed:

  Assumed survival probability over the step.

- benefit:

  Benefit paid at start of step. Default 0.

- h:

  Step length. Default 1.

## Value

Numeric vector.

## Examples

``` r
GI_cont(Vt = 10, Vt1 = 11, P = 1, delta_actual = 0.05, p_assumed = 0.99)
#> [1] 0.6739821
```

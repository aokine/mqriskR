# Interest gain for a continuous-style recursion

Evaluates the one-step gain using the actual force of interest and the
assumed survival probability.

## Usage

``` r
GI_cont(Vt, Vt1, P, delta_actual, p_assumed, benefit = 0, h = 1)
```

## Arguments

- Vt:

  Reserve at time \`t\`.

- Vt1:

  Reserve at time \`t + h\`.

- P:

  Premium rate.

- delta_actual:

  Actual force of interest.

- p_assumed:

  Assumed survival probability over the step.

- benefit:

  Benefit paid at the start of the step.

- h:

  Positive step length.

## Value

Numeric vector of interest gain values.

## Examples

``` r
GI_cont(
  Vt = 10,
  Vt1 = 11,
  P = 1,
  delta_actual = 0.05,
  p_assumed = 0.99
)
#> [1] 0.6739821
```

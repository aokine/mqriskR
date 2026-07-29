# Total gain for a continuous-style one-step recursion

Computes the amount accumulated during a step, less the expected reserve
required at the end of the step.

## Usage

``` r
GT_cont(Vt, Vt1, P, delta_actual, p_actual, benefit = 0, h = 1)
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

- p_actual:

  Actual survival probability over the step.

- benefit:

  Benefit paid at the start of the step.

- h:

  Positive step length.

## Value

Numeric vector of gain values.

## Details

Reserves, premium rates, benefits, and forces of interest may be
negative when such values are meaningful for the application. The step
length must be positive, and the survival probability must lie in
\\\[0,1\]\\.

## Examples

``` r
GT_cont(
  Vt = 10,
  Vt1 = 11,
  P = 1,
  delta_actual = 0.05,
  p_actual = 0.99
)
#> [1] 0.6739821
```

# Gain or loss in a multiple-decrement model

Computes the gain or loss expression from Section 14.6.

## Usage

``` r
gain_loss_md(
  Vt,
  G,
  r,
  e,
  i,
  b1,
  b2,
  s1 = 0,
  s2 = 0,
  q1,
  q2,
  Vt1,
  year_end_cause2 = FALSE,
  q1prime = NULL,
  q2prime = NULL
)
```

## Arguments

- Vt:

  Gross premium reserve at time \\t\\.

- G:

  Gross premium for the year.

- r:

  Percent-of-premium expense factor.

- e:

  Fixed expense at the beginning of the year.

- i:

  Earned interest rate.

- b1:

  Cause 1 benefit.

- b2:

  Cause 2 benefit.

- s1:

  Claim settlement expense for Cause 1.

- s2:

  Claim settlement expense for Cause 2.

- q1:

  Cause 1 decrement probability.

- q2:

  Cause 2 decrement probability.

- Vt1:

  Gross premium reserve at time \\t+1\\.

- year_end_cause2:

  Logical; if `TRUE`, use the year-end Cause 2 form.

- q1prime:

  Single-decrement Cause 1 probability for the year-end Cause 2 case.

- q2prime:

  Single-decrement Cause 2 probability for the year-end Cause 2 case.

## Value

A numeric scalar.

## Details

With within-year decrement probabilities, the function evaluates \$\$
\[{}\_{t}V^G + G(1-r) - e\](1+i) - \left\[(b^{(1)}+s^{(1)})q^{(1)} +
(b^{(2)}+s^{(2)})q^{(2)} + p^{(\tau)} {}\_{t+1}V^G \right\] \$\$

If `year_end_cause2 = TRUE`, the Cause 2 decrement is treated as
occurring only at year end, matching Equation (14.30).

## Examples

``` r
gain_loss_md(
  Vt = 115.00, G = 16, r = 0, e = 3, i = 0.06,
  b1 = 1000, b2 = 110, s1 = 0, s2 = 0,
  q1 = 0.01, q2 = 0.10, Vt1 = 128.83
)
#> [1] 0.0213
```

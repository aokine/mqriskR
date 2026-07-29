# Gain or loss in a two-cause multiple-decrement model

Computes one-year gain or loss under simultaneous within-year decrement
probabilities or an ordered case in which Cause 2 occurs at year-end.
Numeric arguments may be scalar or compatible vectors.

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

  Gross reserve at the beginning of the year.

- G:

  Gross premium.

- r:

  Percent-of-premium expense rate.

- e:

  Fixed beginning-of-year expense.

- i:

  Earned effective annual interest rate.

- b1:

  Cause 1 benefit.

- b2:

  Cause 2 benefit.

- s1:

  Cause 1 settlement expense.

- s2:

  Cause 2 settlement expense.

- q1:

  Cause 1 decrement probability.

- q2:

  Cause 2 decrement probability.

- Vt1:

  Gross reserve at the end of the year.

- year_end_cause2:

  Whether Cause 2 occurs only at year-end.

- q1prime:

  Single-decrement Cause 1 probability for the ordered case.

- q2prime:

  Single-decrement Cause 2 probability for the ordered case.

## Value

A numeric vector of gains or losses.

## Examples

``` r
gain_loss_md(
  Vt = 115, G = 16, r = 0, e = 3, i = 0.06,
  b1 = 1000, b2 = 110, q1 = 0.01, q2 = 0.10, Vt1 = 128.83
)
#> [1] 0.0213
```

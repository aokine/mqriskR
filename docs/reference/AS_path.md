# Projected asset-share path for two decrement causes

Convenience interface to
[`AS_path_md()`](https://aokine.github.io/mqriskR/reference/AS_path_md.md)
for two causes.

## Usage

``` r
AS_path(AS0, G, r, e, b1, b2, q1, q2, p_tau, i, b3 = NULL)
```

## Arguments

- AS0:

  Initial asset share.

- G:

  Premium amount by policy year.

- r:

  Percent-of-premium expense rate by policy year.

- e:

  Fixed expense by policy year.

- b1:

  Cause 1 benefit by policy year.

- b2:

  Cause 2 benefit by policy year.

- q1:

  Cause 1 probability by policy year.

- q2:

  Cause 2 probability by policy year.

- p_tau:

  In-force probability by policy year.

- i:

  Effective annual interest rate by policy year.

- b3:

  Survival benefit by policy year.

## Value

A data frame with policy year `k` and asset share `AS`.

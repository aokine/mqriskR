# Credited rates from index growth rates

Applies a participation rate, floor, cap, and optional margin to raw
index growth rates.

## Usage

``` r
i_credit_eiul(
  i_raw,
  part = 1,
  floor = 0,
  cap = Inf,
  margin = 0,
  margin_after_participation = TRUE
)
```

## Arguments

- i_raw:

  Numeric vector of raw index growth rates.

- part:

  Nonnegative scalar participation rate.

- floor:

  Scalar minimum credited rate.

- cap:

  Scalar maximum credited rate. The default is `Inf`.

- margin:

  Nonnegative scalar index margin.

- margin_after_participation:

  Logical scalar. If `TRUE`, the margin is subtracted after applying
  participation; otherwise it is subtracted before applying
  participation.

## Value

A numeric vector of credited rates.

## Examples

``` r
raw <- iP_eiul(c(1000, 1050, 1200, 1100))
i_credit_eiul(raw, part = 1.10, floor = 0.01, cap = 0.10)
#> [1] 0.055 0.100 0.010
i_credit_eiul(raw)
#> [1] 0.0500000 0.1428571 0.0000000
```

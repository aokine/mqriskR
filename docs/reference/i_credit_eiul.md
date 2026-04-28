# Credited rates from raw index growth rates

Applies participation, floor, cap, and optional margin to raw index
growth rates for an indexed universal life contract.

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

  Participation rate.

- floor:

  Minimum credited rate.

- cap:

  Maximum credited rate.

- margin:

  Index margin. Defaults to 0.

- margin_after_participation:

  Logical; if `TRUE`, subtract the margin after applying the
  participation rate.

## Value

Numeric vector of credited rates.

## Examples

``` r
raw <- iP_eiul(c(1000, 1050, 1200, 1100))
i_credit_eiul(raw, part = 1.10, floor = 0.01, cap = 0.10)
#> [1] 0.055 0.100 0.010
```

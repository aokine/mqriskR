# Zeroized reserves for a discrete death-benefit contract

Computes reserves backward by setting expected profit in each policy
year equal to zero. Negative reserves may optionally be floored at zero.

## Usage

``` r
V_zeroized(qx, i, G, benefit, r = 0, e = 0, V_terminal = 0, floor_zero = TRUE)
```

## Arguments

- qx:

  Mortality probability by policy year.

- i:

  Annual effective interest rate by policy year. Values must be greater
  than `-1`.

- G:

  Gross premium by policy year.

- benefit:

  Death benefit by policy year.

- r:

  Percent-of-premium expense rate by policy year. Values must lie in
  `[0, 1]`.

- e:

  Fixed expense by policy year.

- V_terminal:

  Nonnegative scalar terminal reserve.

- floor_zero:

  Logical scalar. If `TRUE`, negative reserves are replaced by zero.

## Value

A named numeric vector of length `length(qx) + 1`.

## Examples

``` r
V_zeroized(
  qx = c(0.015, 0.017, 0.019, 0.021, 0.024),
  i = 0.06,
  G = 19279,
  benefit = 1000000,
  e = 240
)
#>       V0       V1       V2       V3       V4       V5 
#>    0.000    0.000 2679.540 4099.544 3602.509    0.000 
```

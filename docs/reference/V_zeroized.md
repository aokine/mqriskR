# Zeroized reserves for a discrete death-only contract

Computes the zeroized reserve sequence by backward recursion, setting
negative reserves equal to zero.

## Usage

``` r
V_zeroized(qx, i, G, benefit, r = 0, e = 0, V_terminal = 0, floor_zero = TRUE)
```

## Arguments

- qx:

  Mortality vector.

- i:

  Interest-rate vector.

- G:

  Gross premium vector.

- benefit:

  Death-benefit vector.

- r:

  Percent-of-premium expense vector.

- e:

  Fixed-expense vector.

- V_terminal:

  Terminal reserve. Defaults to 0.

- floor_zero:

  Logical; if `TRUE`, negative reserves are reset to 0.

## Value

Numeric vector of zeroized reserves of length \\n+1\\.

## Details

For a death-only contract with no settlement expense and no second
decrement, the recursion sets \$\$ Pr\_{t+1} = ({}\_tV^Z +
G\_{t+1}(1-r\_{t+1}) - e\_{t+1})(1+i\_{t+1}) - \[Bq\_{x+t} +
{}\_{t+1}V^Z p\_{x+t}\] \$\$ equal to zero, solving backward for
\\{}\_tV^Z\\.

## Examples

``` r
V_zeroized(
  qx = c(.015, .017, .019, .021, .024),
  i = 0.06,
  G = 19279,
  benefit = 1000000,
  e = 240
)
#>       V0       V1       V2       V3       V4       V5 
#>    0.000    0.000 2679.540 4099.544 3602.509    0.000 
```

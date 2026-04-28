# Reserve derivatives for the disability model with recovery

Computes the right-hand sides of the coupled Thiele differential
equations in Equations (14.25) and (14.26) for the healthy-life reserve
\\{}\_{t}\overline{V}^{(0)}\\ and the disabled-life reserve
\\{}\_{t}\overline{V}^{(1)}\\.

## Usage

``` r
thiele_dVdt_01(t, V0, V1, delta, Pbar, B, R, mu01, mu02, mu10, mu12)
```

## Arguments

- t:

  Time.

- V0:

  Value of \\{}\_{t}\overline{V}^{(0)}\\.

- V1:

  Value of \\{}\_{t}\overline{V}^{(1)}\\.

- delta:

  Force of interest.

- Pbar:

  Continuous premium rate.

- B:

  Death benefit.

- R:

  Continuous disability income rate.

- mu01:

  Function of time returning \\\mu\_{x+t}^{01}\\.

- mu02:

  Function of time returning \\\mu\_{x+t}^{02}\\.

- mu10:

  Function of time returning \\\mu\_{x+t}^{10}\\.

- mu12:

  Function of time returning \\\mu\_{x+t}^{12}\\.

## Value

A named numeric vector with components `dV0` and `dV1`.

## Details

The equations are \$\$ \frac{d}{dt}{}\_{t}\overline{V}^{(0)} =
\overline{P} + \delta {}\_{t}\overline{V}^{(0)} -
\mu\_{x+t}^{02}(B-{}\_{t}\overline{V}^{(0)}) -
\mu\_{x+t}^{01}({}\_{t}\overline{V}^{(1)}-{}\_{t}\overline{V}^{(0)})
\$\$ and \$\$ \frac{d}{dt}{}\_{t}\overline{V}^{(1)} = \delta
{}\_{t}\overline{V}^{(1)} - R -
\mu\_{x+t}^{12}(B-{}\_{t}\overline{V}^{(1)}) -
\mu\_{x+t}^{10}({}\_{t}\overline{V}^{(0)}-{}\_{t}\overline{V}^{(1)})
\$\$

## Examples

``` r
mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

thiele_dVdt_01(
  t = 2.0, V0 = 0, V1 = 0,
  delta = 0.04, Pbar = 446.95,
  B = 1000, R = 1000,
  mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
)
#>      dV0      dV1 
#>   246.95 -1450.00 
```

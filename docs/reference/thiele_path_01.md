# Backward reserve path for the disability model with recovery

Computes the backward Euler reserve path for the healthy-life reserve
\\{}\_{t}\overline{V}^{(0)}\\ and disabled-life reserve
\\{}\_{t}\overline{V}^{(1)}\\ using Equations (14.27) and (14.28).

## Usage

``` r
thiele_path_01(
  h,
  n,
  delta,
  Pbar,
  B,
  R,
  mu01,
  mu02,
  mu10,
  mu12,
  V0_n = 0,
  V1_n = 0
)
```

## Arguments

- h:

  Step size.

- n:

  Final time.

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

- V0_n:

  Terminal value of \\{}\_{n}\overline{V}^{(0)}\\.

- V1_n:

  Terminal value of \\{}\_{n}\overline{V}^{(1)}\\.

## Value

A data frame with columns `t`, `tV0`, and `tV1`.

## Examples

``` r
mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

thiele_path_01(
  h = 0.10, n = 2.0, delta = 0.04, Pbar = 446.95,
  B = 1000, R = 1000,
  mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
)
#>      t        tV0       tV1
#> 1  0.0  -5.386519 1283.8209
#> 2  0.1  -7.435183 1257.5752
#> 3  0.2 -10.264825 1229.6573
#> 4  0.3 -13.816667 1199.8895
#> 5  0.4 -18.020269 1168.0679
#> 6  0.5 -22.790964 1133.9575
#> 7  0.6 -28.026787 1097.2878
#> 8  0.7 -33.604789 1057.7466
#> 9  0.8 -39.376622 1014.9732
#> 10 0.9 -45.163231  968.5502
#> 11 1.0 -50.748473  917.9942
#> 12 1.1 -55.871450  862.7435
#> 13 1.2 -60.217265  802.1453
#> 14 1.3 -63.405885  735.4389
#> 15 1.4 -64.978700  661.7368
#> 16 1.5 -64.382289  580.0010
#> 17 1.6 -60.948789  489.0159
#> 18 1.7 -53.872127  387.3551
#> 19 1.8 -42.179215  273.3415
#> 20 1.9 -24.695000  145.0000
#> 21 2.0   0.000000    0.0000
```

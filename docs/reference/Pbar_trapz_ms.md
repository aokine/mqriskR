# Continuous premium approximation \\\overline{P}\\ by trapezoidal rule

Approximates the annual continuous premium in the disability model
allowing for recovery, as in Example 14.18.

## Usage

``` r
Pbar_trapz_ms(t, tp00, tp01, delta, mu02, mu12, B02 = 1, B12 = 1, R = 0)
```

## Arguments

- t:

  Numeric vector of time points.

- tp00:

  Numeric vector of values \\{}\_{t}p\_{x}^{00}\\.

- tp01:

  Numeric vector of values \\{}\_{t}p\_{x}^{01}\\.

- delta:

  Force of interest.

- mu02:

  Function of time returning \\\mu\_{x+t}^{02}\\.

- mu12:

  Function of time returning \\\mu\_{x+t}^{12}\\.

- B02:

  Benefit payable on death while healthy.

- B12:

  Benefit payable on death while disabled.

- R:

  Continuous income rate while disabled.

## Value

A numeric scalar.

## Details

The numerator is \$\$ \int v^t
\left\[{}\_{t}p\_{x}^{00}\mu\_{x+t}^{02}B^{02} +
{}\_{t}p\_{x}^{01}\mu\_{x+t}^{12}B^{12} + {}\_{t}p\_{x}^{01}R \right\]
dt \$\$ and the denominator is \$\$ \int v^t {}\_{t}p\_{x}^{00} dt \$\$

## Examples

``` r
mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

ex1410 <- tp00_tp01_euler(
  h = 0.10, n = 2.0,
  mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
)

Pbar_trapz_ms(
  t = ex1410$t,
  tp00 = ex1410$tp00,
  tp01 = ex1410$tp01,
  delta = 0.04,
  mu02 = mu02,
  mu12 = mu12,
  B02 = 1000,
  B12 = 1000,
  R = 1000
)
#> [1] 446.9451
```

# Continuous premium approximation in a disability model

Approximates a continuous premium rate by trapezoidal integration.

## Usage

``` r
Pbar_trapz_ms(t, tp00, tp01, delta, mu02, mu12, B02 = 1, B12 = 1, R = 0)
```

## Arguments

- t:

  Strictly increasing nonnegative time points.

- tp00:

  Healthy-state probabilities.

- tp01:

  Disabled-state probabilities.

- delta:

  Force of interest.

- mu02:

  Healthy-to-deceased intensity function.

- mu12:

  Disabled-to-deceased intensity function.

- B02:

  Benefit on death while healthy.

- B12:

  Benefit on death while disabled.

- R:

  Continuous disability income rate.

## Value

A numeric scalar.

## Examples

``` r
mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

probs <- tp00_tp01_euler(
  h = 0.10,
  n = 2,
  mu01 = mu01,
  mu02 = mu02,
  mu10 = mu10,
  mu12 = mu12
)

Pbar_trapz_ms(
  t = probs$t,
  tp00 = probs$tp00,
  tp01 = probs$tp01,
  delta = 0.04,
  mu02 = mu02,
  mu12 = mu12,
  B02 = 1000,
  B12 = 1000,
  R = 1000
)
#> [1] 446.9451
```

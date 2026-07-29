# Backward reserve path for a disability model with recovery

Computes a backward Euler reserve path from terminal healthy-state and
disabled-state reserves.

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

  Positive step size.

- n:

  Nonnegative final time.

- delta:

  Force of interest.

- Pbar:

  Continuous premium rate.

- B:

  Death benefit.

- R:

  Continuous disability income rate.

- mu01:

  Healthy-to-disabled intensity function.

- mu02:

  Healthy-to-deceased intensity function.

- mu10:

  Disabled-to-healthy intensity function.

- mu12:

  Disabled-to-deceased intensity function.

- V0_n:

  Terminal healthy-state reserve.

- V1_n:

  Terminal disabled-state reserve.

## Value

A data frame with columns `t`, `tV0`, and `tV1`.

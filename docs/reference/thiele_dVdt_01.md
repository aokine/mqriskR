# Reserve derivatives for a disability model with recovery

Computes the coupled Thiele reserve derivatives. Numeric arguments may
be scalar or compatible vectors.

## Usage

``` r
thiele_dVdt_01(t, V0, V1, delta, Pbar, B, R, mu01, mu02, mu10, mu12)
```

## Arguments

- t:

  Time.

- V0:

  Healthy-state reserve.

- V1:

  Disabled-state reserve.

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

## Value

A named vector for scalar input or a two-column matrix for vectorized
input.

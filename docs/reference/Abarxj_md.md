# Continuous multiple-decrement insurance present value

Approximates the actuarial present value of a benefit payable at the
moment of decrement using the trapezoidal rule.

## Usage

``` r
Abarxj_md(t, ptau, muj, delta, benefit = 1)
```

## Arguments

- t:

  Strictly increasing nonnegative time points.

- ptau:

  In-force probabilities at the supplied time points.

- muj:

  Cause-specific decrement intensities at the supplied time points.

- delta:

  Force of interest. May be scalar or vector.

- benefit:

  Benefit payable on decrement. May be scalar or vector.

## Value

A numeric vector of actuarial present values.

## Examples

``` r
t <- seq(0, 20, by = 0.1)
ptau <- exp(-0.012 * t)
muj <- rep(0.002, length(t))
Abarxj_md(t, ptau, muj, delta = 0.05, benefit = 2000)
#> [1] 45.84633
```

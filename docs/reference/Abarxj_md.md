# Continuous multiple-decrement insurance APV \\\overline{A}\_{x}^{(j)}\\

Computes the actuarial present value of a benefit payable at the moment
of decrement by Cause \\j\\, matching Equation (14.4) in Chapter 14.

## Usage

``` r
Abarxj_md(t, ptau, muj, delta, benefit = 1)
```

## Arguments

- t:

  Numeric vector of time points.

- ptau:

  Numeric vector of values \\{}\_{t}p\_{x}^{(\tau)}\\.

- muj:

  Numeric vector of values \\\mu\_{x+t}^{(j)}\\.

- delta:

  Force of interest.

- benefit:

  Benefit amount payable on decrement by Cause \\j\\.

## Value

A numeric scalar.

## Details

The integral is evaluated numerically by the trapezoidal rule: \$\$
\overline{A}\_{x}^{(j)} = \int_0^T v^t {}\_{t}p\_{x}^{(\tau)}
\mu\_{x+t}^{(j)} dt \$\$

## Examples

``` r
t <- seq(0, 20, by = 0.01)
ptau <- exp(-0.012 * t)
mu_ac <- rep(0.002, length(t))
Abarxj_md(t, ptau, mu_ac, delta = 0.05, benefit = 2000)
#> [1] 45.84618
```

# Discrete multiple-decrement insurance present value

Computes the actuarial present value of a benefit payable at the end of
the year of decrement from a specified cause.

## Usage

``` r
Axj_md(qj, ptau, i, benefit = 1)
```

## Arguments

- qj:

  Cause-specific decrement probabilities by policy year.

- ptau:

  In-force probabilities at the beginning of each policy year.

- i:

  Effective annual interest rate. May be scalar or vector.

- benefit:

  Benefit payable on decrement. May be scalar or vector.

## Value

A numeric vector of actuarial present values.

## Examples

``` r
qj <- c(0.02, 0.02, 0.02, 0.02, 0.02)
ptau <- c(1, 0.95, 0.89, 0.82, 0.74)
Axj_md(qj, ptau, i = 0.06, benefit = 1000)
#> [1] 74.77284
```

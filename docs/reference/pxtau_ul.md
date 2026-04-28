# One-year persistency rates for universal life

Computes the one-year persistency rates \\p_x^{(\tau)}\\ under either
Equation (16.14) or Equation (16.15).

## Usage

``` r
pxtau_ul(qd, qw, year_end_withdrawal = TRUE)
```

## Arguments

- qd:

  Mortality probabilities.

- qw:

  Withdrawal probabilities.

- year_end_withdrawal:

  Logical; if `TRUE`, use \\(1-q^{(d)})(1-q^{(w)})\\. Otherwise use
  \\1-q^{(d)}-q^{(w)}\\.

## Value

Numeric vector.

## Examples

``` r
qd <- c(.001, .002, .003)
qw <- c(.02, .02, .03)
pxtau_ul(qd, qw)
#> [1] 0.97902 0.97804 0.96709
```

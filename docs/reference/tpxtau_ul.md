# Cumulative persistency to the end of each policy year

Computes \\{}\_tp_x^{(\tau)}\\ from the one-year persistency rates.

## Usage

``` r
tpxtau_ul(qd, qw, year_end_withdrawal = TRUE)
```

## Arguments

- qd:

  Mortality probabilities.

- qw:

  Withdrawal probabilities.

- year_end_withdrawal:

  Logical; if `TRUE`, use Equation (16.15).

## Value

Numeric vector.

## Examples

``` r
qd <- c(.001, .002, .003)
qw <- c(.02, .02, .03)
tpxtau_ul(qd, qw)
#> [1] 0.9790200 0.9575207 0.9260087
```

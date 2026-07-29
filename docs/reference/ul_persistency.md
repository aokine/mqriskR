# Universal life persistency probabilities

`pxtau_ul()` computes one-year persistency probabilities under mortality
and withdrawal.

## Usage

``` r
pxtau_ul(qd, qw, year_end_withdrawal = TRUE)

tpxtau_ul(qd, qw, year_end_withdrawal = TRUE)
```

## Arguments

- qd:

  Mortality probabilities.

- qw:

  Withdrawal probabilities.

- year_end_withdrawal:

  Logical scalar. If `TRUE`, withdrawal is modeled at year-end and
  persistency is \\(1-q^{(d)})(1-q^{(w)})\\. Otherwise persistency is
  \\1-q^{(d)}-q^{(w)}\\.

## Value

A numeric vector.

## Details

`tpxtau_ul()` computes cumulative persistency through the end of each
policy year.

## Examples

``` r
qd <- c(0.001, 0.002, 0.003)
qw <- c(0.02, 0.02, 0.03)

pxtau_ul(qd, qw)
#> [1] 0.97902 0.97804 0.96709
tpxtau_ul(qd, qw)
#> [1] 0.9790200 0.9575207 0.9260087
```

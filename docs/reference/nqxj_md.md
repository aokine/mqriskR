# Cause-specific multiple-decrement probability from a table

Computes the probability of decrement from cause `j` within `n` years:
\$\$ {}\_nq_x^{(j)} = \sum\_{k=0}^{n-1} {}\_kp_x^{(\tau)}q\_{x+k}^{(j)}.
\$\$

## Usage

``` r
nqxj_md(tbl, x, n, j)
```

## Arguments

- tbl:

  A multiple-decrement table produced by
  [`md_table()`](https://aokine.github.io/mqriskR/reference/md_table.md).

- x:

  Starting age or duration. May be scalar or vector.

- n:

  Nonnegative integer term. May be scalar or vector.

- j:

  Positive integer cause index. May be scalar or vector.

## Value

A numeric vector of cause-specific decrement probabilities.

## Examples

``` r
ages <- 45:50
qmat <- cbind(
  q1 = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
  q2 = rep(0.100, 6)
)
tbl <- md_table(ages, qmat, radix = 1000)

nqxj_md(tbl, x = 46, n = 2, j = 1)
#> [1] 0.023544
```

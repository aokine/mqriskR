# Multiple-decrement survival probability from a table

Computes \$\${}\_np_x^{(\tau)} =\prod\_{k=0}^{n-1}p\_{x+k}^{(\tau)}.\$\$

## Usage

``` r
npxtau_md(tbl, x, n)
```

## Arguments

- tbl:

  A multiple-decrement table produced by
  [`md_table()`](https://aokine.github.io/mqriskR/reference/md_table.md).

- x:

  Starting age or duration. May be scalar or vector.

- n:

  Nonnegative integer term. May be scalar or vector.

## Value

A numeric vector of survival probabilities.

## Examples

``` r
ages <- 45:50
qmat <- cbind(
  q1 = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
  q2 = rep(0.100, 6)
)
tbl <- md_table(ages, qmat, radix = 1000)

npxtau_md(tbl, x = 46, n = 3)
#> [1] 0.6978632
```

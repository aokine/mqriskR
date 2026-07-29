# Total multiple-decrement probability from a table

Computes \$\${}\_nq_x^{(\tau)}=1-{}\_np_x^{(\tau)}.\$\$

## Usage

``` r
nqxtau_md(tbl, x, n)
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

A numeric vector of total decrement probabilities.

## Examples

``` r
ages <- 45:50
qmat <- cbind(
  q1 = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
  q2 = rep(0.100, 6)
)
tbl <- md_table(ages, qmat, radix = 1000)

nqxtau_md(tbl, x = 46, n = 2)
#> [1] 0.212344
```

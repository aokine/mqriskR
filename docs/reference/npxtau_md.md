# \\{}\_n p_x^{(\tau)}\\ from a multiple-decrement table

\\{}\_n p_x^{(\tau)}\\ from a multiple-decrement table

## Usage

``` r
npxtau_md(tbl, x, n)
```

## Arguments

- tbl:

  Output from md_table().

- x:

  Starting age.

- n:

  Number of years.

## Value

Numeric scalar.

## Examples

``` r
x <- 45:50
qmat <- cbind(q1 = c(.011, .012, .013, .014, .015, .016), q2 = rep(.1, 6))
tbl <- md_table(x, qmat, radix = 1000)
npxtau_md(tbl, x = 46, n = 3)
#> [1] 0.6978632
```

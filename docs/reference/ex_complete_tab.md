# Complete expectation of life from a life table

Computes \\\overset{\circ}{e}\_x = \int_0^\infty {}\_t p_x \\ dt\\ using
a within-year assumption: `"udd"`, `"cf"`, or `"balducci"`.

## Usage

``` r
ex_complete_tab(tbl, x, assumption = c("udd", "cf", "balducci"))
```

## Arguments

- tbl:

  A life_table object.

- x:

  Numeric vector of integer ages.

- assumption:

  One of `"udd"`, `"cf"`, `"balducci"`.

## Value

Numeric vector of complete expectations \\\overset{\circ}{e}\_x\\.

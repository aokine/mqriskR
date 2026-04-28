# Temporary complete expectation of life from a life table

Computes \\\overset{\circ}{e}\_{x:\overline{n}\|} = \int_0^n {}\_t p_x
\\ dt\\ using a within-year assumption: `"udd"`, `"cf"`, or
`"balducci"`.

## Usage

``` r
ex_temp_complete_tab(tbl, x, n, assumption = c("udd", "cf", "balducci"))
```

## Arguments

- tbl:

  A life_table object.

- x:

  Numeric vector of integer ages.

- n:

  Numeric vector of nonnegative numbers.

- assumption:

  One of `"udd"`, `"cf"`, `"balducci"`.

## Value

Numeric vector of temporary complete expectations.

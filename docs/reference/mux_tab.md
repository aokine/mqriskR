# Fractional force of mortality from a life table

Computes \\\mu\_{x+t}\\ under UDD, constant force, or Balducci.

## Usage

``` r
mux_tab(tbl, x, t, assumption = c("udd", "cf", "balducci"))
```

## Arguments

- tbl:

  A life_table object.

- x:

  Numeric vector of integer ages.

- t:

  Numeric vector of fractional durations with \\0 \le t \le 1\\.

- assumption:

  One of `"udd"`, `"cf"`, `"balducci"`.

## Value

Numeric vector of \\\mu\_{x+t}\\ values.

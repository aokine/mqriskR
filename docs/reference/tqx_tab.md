# Fractional failure probability from a life table

Computes \\{}\_t q_x = 1 - {}\_t p_x\\ for \\0 \le t \le 1\\.

## Usage

``` r
tqx_tab(tbl, x, t, assumption = c("udd", "cf", "balducci"))
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

Numeric vector of \\{}\_t q_x\\ values.

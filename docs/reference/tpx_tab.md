# Fractional survival probability from a life table

Computes \\{}\_t p_x\\ for \\0 \le t \le 1\\ from a discrete life table
under one of the standard Chapter 6 assumptions: UDD, constant force, or
Balducci.

## Usage

``` r
tpx_tab(tbl, x, t, assumption = c("udd", "cf", "balducci"))
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

Numeric vector of \\{}\_t p_x\\ values.

# Fractional conditional density from a life table

Computes the conditional density \\f_x(t \mid T_0 \> x) = {}\_t p_x
\mu\_{x+t}\\ for \\0 \< t \< 1\\.

## Usage

``` r
fx_tab(tbl, x, t, assumption = c("udd", "cf", "balducci"))
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

Numeric vector of conditional density values.

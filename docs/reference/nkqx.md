# Curtate death probability from a life table

Computes \\{}\_{k\|} q_x = {}\_k p_x - {}\_{k+1} p_x\\

## Usage

``` r
nkqx(tbl, x, k)
```

## Arguments

- tbl:

  A life_table object.

- x:

  Numeric vector of ages.

- k:

  Nonnegative integer.

## Value

Numeric vector of \\{}\_{k\|} q_x\\ values.

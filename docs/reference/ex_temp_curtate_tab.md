# Temporary curtate expectation of life from a life table

Computes \\e\_{x:\overline{n}\|} = \sum\_{k=1}^{n} {}\_k p_x\\ for
integer n in the discrete tabular setting.

## Usage

``` r
ex_temp_curtate_tab(tbl, x, n)
```

## Arguments

- tbl:

  A life_table object.

- x:

  Numeric vector of integer ages.

- n:

  Numeric vector of nonnegative integers.

## Value

Numeric vector of temporary curtate expectations.

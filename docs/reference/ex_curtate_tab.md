# Curtate expectation of life from a life table

Computes the curtate expectation of life \\e_x = \sum\_{k=1}^\infty
{}\_k p_x\\ in the discrete tabular setting.

## Usage

``` r
ex_curtate_tab(tbl, x)
```

## Arguments

- tbl:

  A life_table object.

- x:

  Numeric vector of integer ages.

## Value

Numeric vector of curtate expectations \\e_x\\.

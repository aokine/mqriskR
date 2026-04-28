# Deferred death probability from a life table

Computes the probability that a life aged x survives n years and then
dies within the following m years: \\{}\_{n\|m} q_x = {}\_n p_x \cdot
{}\_m q\_{x+n}\\

## Usage

``` r
nmxq(tbl, x, n, m)
```

## Arguments

- tbl:

  A life_table object.

- x:

  Numeric vector of ages.

- n:

  Nonnegative integer deferred period.

- m:

  Nonnegative integer subsequent period.

## Value

Numeric vector of \\{}\_{n\|m} q_x\\ values.

## Details

This function is for integer n and m in the discrete tabular setting.

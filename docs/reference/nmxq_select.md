# Deferred select-life death probability

Computes \\{}\_{n\|m} q\_{\[x\]+t} = {}\_n p\_{\[x\]+t} \cdot {}\_m
q\_{\[x\]+t+n}\\.

## Usage

``` r
nmxq_select(tbl, x_sel, t, n, m)
```

## Arguments

- tbl:

  A select_life_table object.

- x_sel:

  Numeric vector of ages at selection.

- t:

  Numeric vector of current durations since selection.

- n:

  Numeric vector of nonnegative integer deferred periods.

- m:

  Numeric vector of nonnegative integer death windows.

## Value

Numeric vector of deferred death probabilities.

# Select-life death probability

Computes \\{}\_n q\_{\[x\]+t} = 1 - {}\_n p\_{\[x\]+t}\\.

## Usage

``` r
nqx_select(tbl, x_sel, t, n)
```

## Arguments

- tbl:

  A select_life_table object.

- x_sel:

  Numeric vector of ages at selection.

- t:

  Numeric vector of current durations since selection.

- n:

  Numeric vector of nonnegative integer future durations.

## Value

Numeric vector of death probabilities.

# Extract select-table survivor value

Returns \\l\_{\[x\]+t}\\ from a select life table.

## Usage

``` r
lx_select(tbl, x_sel, t)
```

## Arguments

- tbl:

  A select_life_table object.

- x_sel:

  Numeric vector of ages at selection.

- t:

  Numeric vector of durations since selection.

## Value

Numeric vector of survivor values.

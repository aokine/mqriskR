# Select-life survival probability

Computes \\{}\_n p\_{\[x\]+t} = l\_{\[x\]+t+n} / l\_{\[x\]+t}\\

## Usage

``` r
npx_select(tbl, x_sel, t, n)
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

Numeric vector of survival probabilities.

## Details

in the discrete select-table setting.

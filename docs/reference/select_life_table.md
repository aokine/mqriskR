# Construct a select life table

Builds a select-life-table object from vectors of selection age,
duration since selection, attained age, and survivor values.

## Usage

``` r
select_life_table(x_sel, duration, attained_age, lx)
```

## Arguments

- x_sel:

  Numeric vector of ages at selection.

- duration:

  Numeric vector of durations since selection.

- attained_age:

  Numeric vector of attained ages.

- lx:

  Numeric vector of select-table survivor values.

## Value

A data.frame with class `"select_life_table"`.

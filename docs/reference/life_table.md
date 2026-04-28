# Construct a life table

Build a discrete life table from one of lx, qx, px, or S0.

## Usage

``` r
life_table(x, lx = NULL, qx = NULL, px = NULL, S0 = NULL, radix = 1e+05)
```

## Arguments

- x:

  Numeric vector of ages.

- lx:

  Numeric vector of l_x values.

- qx:

  Numeric vector of q_x values.

- px:

  Numeric vector of p_x values.

- S0:

  Numeric vector of S_0(x) values.

- radix:

  Radix used when converting S0 to lx, or when building from qx/px.

## Value

A data.frame with class "life_table".

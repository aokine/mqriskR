# Construct life-table values from q_x values

Builds life-table survivor values recursively from \\l\_{x+1} = l_x
(1-q_x)\\, starting from a chosen radix.

## Usage

``` r
qx_to_lx(qx, radix = 1e+05)
```

## Arguments

- qx:

  Numeric vector of one-year death probabilities \\q_x\\.

- radix:

  Positive radix \\l_0\\.

## Value

Numeric vector of \\l_x\\ values of length `length(qx)+1`.

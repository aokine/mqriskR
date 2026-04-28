# Construct life-table values from p_x values

Builds life-table survivor values recursively from \\l\_{x+1} = l_x
p_x\\, starting from a chosen radix.

## Usage

``` r
px_to_lx(px, radix = 1e+05)
```

## Arguments

- px:

  Numeric vector of one-year survival probabilities \\p_x\\.

- radix:

  Positive radix \\l_0\\.

## Value

Numeric vector of \\l_x\\ values of length `length(px)+1`.

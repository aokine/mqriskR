# Construct a multiple-decrement table

Constructs a discrete multiple-decrement table from cause-specific
one-year decrement probabilities.

## Usage

``` r
md_table(x, qxj, radix = 1e+05)
```

## Arguments

- x:

  Integer vector of consecutive ages or durations.

- qxj:

  Numeric matrix or data frame of cause-specific decrement
  probabilities. Rows correspond to values in `x`, and columns
  correspond to causes.

- radix:

  Initial value of \\l_x^{(\tau)}\\.

## Value

An object of classes `"md_table"` and `"data.frame"` containing age or
duration, cause-specific decrement probabilities, total decrement and
survival probabilities, numbers alive, and cause-specific and total
decrements.

## Examples

``` r
ages <- 45:50
qmat <- cbind(
  withdrawal = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
  retirement = rep(0.100, 6)
)

md_table(ages, qmat, radix = 1000)
#>    x    q1  q2  qtau  ptau      ltau        d1        d2      dtau
#> 1 45 0.011 0.1 0.111 0.889 1000.0000 11.000000 100.00000 111.00000
#> 2 46 0.012 0.1 0.112 0.888  889.0000 10.668000  88.90000  99.56800
#> 3 47 0.013 0.1 0.113 0.887  789.4320 10.262616  78.94320  89.20582
#> 4 48 0.014 0.1 0.114 0.886  700.2262  9.803167  70.02262  79.82578
#> 5 49 0.015 0.1 0.115 0.885  620.4004  9.306006  62.04004  71.34605
#> 6 50 0.016 0.1 0.116 0.884  549.0544  8.784870  54.90544  63.69030
```

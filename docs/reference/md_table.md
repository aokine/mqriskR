# Build a multiple-decrement table

Build a multiple-decrement table

## Usage

``` r
md_table(x, qxj, radix = 1e+05)
```

## Arguments

- x:

  Integer vector of ages or durations.

- qxj:

  Matrix/data.frame of cause-specific decrement probabilities. Rows
  correspond to ages in x, columns correspond to causes.

- radix:

  Starting \\l_x^{(\tau)}\\.

## Value

Data frame containing \\q^{(j)}\\, \\q^{(\tau)}\\, \\p^{(\tau)}\\,
\\l^{(\tau)}\\, \\d^{(j)}\\, and \\d^{(\tau)}\\.

## Examples

``` r
x <- 45:50
qmat <- cbind(
  q1 = c(.011, .012, .013, .014, .015, .016),
  q2 = c(.100, .100, .100, .100, .100, .100)
)
md_table(x, qmat, radix = 1000)
#>    x    q1  q2  qtau  ptau      ltau        d1        d2      dtau
#> 1 45 0.011 0.1 0.111 0.889 1000.0000 11.000000 100.00000 111.00000
#> 2 46 0.012 0.1 0.112 0.888  889.0000 10.668000  88.90000  99.56800
#> 3 47 0.013 0.1 0.113 0.887  789.4320 10.262616  78.94320  89.20582
#> 4 48 0.014 0.1 0.114 0.886  700.2262  9.803167  70.02262  79.82578
#> 5 49 0.015 0.1 0.115 0.885  620.4004  9.306006  62.04004  71.34605
#> 6 50 0.016 0.1 0.116 0.884  549.0544  8.784870  54.90544  63.69030
```

# Account-value path for Type B universal life

Computes the year-by-year account value roll-forward for Type B
universal life using Equations (16.4) and (16.5a).

## Usage

``` r
AV_path_ul_typeB(G, r, e, qx, ic, B, iq = ic, AV0 = 0)
```

## Arguments

- G:

  Premium vector \\G_t\\.

- r:

  Percent-of-premium expense vector \\r_t\\.

- e:

  Fixed expense vector \\e_t\\.

- qx:

  Mortality vector \\q\_{x+t-1}\\.

- ic:

  Credited interest rate vector \\i^c\\.

- B:

  Face amount.

- iq:

  Interest rate vector \\i^q\\ used in cost of insurance. Defaults to
  `ic`.

- AV0:

  Initial account value. Defaults to 0.

## Value

A data frame with columns `t`, `premium`, `net_contribution`, `COI`, and
`AV`.

## Examples

``` r
qx <- c(.00076, .00081, .00085, .00090, .00095)
r <- c(.75, .10, .10, .10, .10)
e <- c(100, 20, 20, 20, 20)
G <- rep(5000, 5)

AV_path_ul_typeB(G = G, r = r, e = e, qx = qx, ic = 0.03, B = 100000)
#>   t premium net_contribution      COI        AV
#> 1 0      NA               NA       NA     0.000
#> 2 1    5000             1150 73.78641  1108.500
#> 3 2    5000             4480 78.64078  5675.155
#> 4 3    5000             4480 82.52427 10374.810
#> 5 4    5000             4480 87.37864 15210.454
#> 6 5    5000             4480 92.23301 20186.168
```

# Type B universal life account-value path

Computes a year-by-year Type B universal life account-value path.
Premiums and expenses are applied at the beginning of each period,
followed by the cost-of-insurance charge and credited interest.

## Usage

``` r
AV_path_ul_typeB(G, r, e, qx, ic, B, iq = ic, AV0 = 0)
```

## Arguments

- G:

  Premium amount by period.

- r:

  Percent-of-premium expense rate by period. Values must lie in
  `[0, 1]`.

- e:

  Fixed expense by period.

- qx:

  Mortality probability by period.

- ic:

  Credited annual effective interest rate by period. Values must be
  greater than `-1`.

- B:

  Face amount by period.

- iq:

  Interest rate used in the cost-of-insurance calculation. Defaults to
  `ic`; values must be greater than `-1`.

- AV0:

  Nonnegative scalar initial account value.

## Value

A data frame containing the policy duration, premium, net contribution,
cost-of-insurance charge, and account value.

## Examples

``` r
qx <- c(0.00076, 0.00081, 0.00085, 0.00090, 0.00095)
r <- c(0.75, rep(0.10, 4))
e <- c(100, rep(20, 4))

AV_path_ul_typeB(
  G = 5000,
  r = r,
  e = e,
  qx = qx,
  ic = 0.03,
  B = 100000
)
#>   t premium net_contribution      COI        AV
#> 1 0      NA               NA       NA     0.000
#> 2 1    5000             1150 73.78641  1108.500
#> 3 2    5000             4480 78.64078  5675.155
#> 4 3    5000             4480 82.52427 10374.810
#> 5 4    5000             4480 87.37864 15210.454
#> 6 5    5000             4480 92.23301 20186.168
```

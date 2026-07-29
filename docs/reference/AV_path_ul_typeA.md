# Type A universal life account-value path

Computes a year-by-year Type A universal life account-value path. The
death benefit is fixed, so the net amount at risk depends on the ending
account value and the roll-forward is solved explicitly each period.

## Usage

``` r
AV_path_ul_typeA(G, r, e, qx, ic, B, iq = ic, AV0 = 0)
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

A data frame containing the policy duration, premium, and account value.

## Examples

``` r
qx <- c(0.00076, 0.00081, 0.00085, 0.00090, 0.00095)
r <- c(0.75, rep(0.10, 4))
e <- c(100, rep(20, 4))

AV_path_ul_typeA(
  G = 5000,
  r = r,
  e = e,
  qx = qx,
  ic = 0.03,
  B = 100000
)
#>   t premium        AV
#> 1 0      NA     0.000
#> 2 1    5000  1109.343
#> 3 2    5000  5680.625
#> 4 3    5000 10389.274
#> 5 4    5000 15239.068
#> 6 5    5000 20234.863
```

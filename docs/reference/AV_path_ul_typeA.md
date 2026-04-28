# Account-value path for Type A universal life

Computes the year-by-year account value roll-forward for Type A
universal life using the explicit form in Equation (16.8), or Equation
(16.9) when \\i^q = i^c\\.

## Usage

``` r
AV_path_ul_typeA(G, r, e, qx, ic, B, iq = ic, AV0 = 0)
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

  Fixed death benefit face amount.

- iq:

  Interest rate vector \\i^q\\ used in cost of insurance. Defaults to
  `ic`.

- AV0:

  Initial account value. Defaults to 0.

## Value

A data frame with columns `t`, `premium`, `AV`.

## Examples

``` r
qx <- c(.00076, .00081, .00085, .00090, .00095)
r <- c(.75, .10, .10, .10, .10)
e <- c(100, 20, 20, 20, 20)
G <- rep(5000, 5)

AV_path_ul_typeA(G = G, r = r, e = e, qx = qx, ic = 0.03, B = 100000)
#>   t premium        AV
#> 1 0      NA     0.000
#> 2 1    5000  1109.343
#> 3 2    5000  5680.625
#> 4 3    5000 10389.274
#> 5 4    5000 15239.068
#> 6 5    5000 20234.863
```

# Spot-rate actuarial present value functions

Chapter 15 functions for actuarial present values when discounting uses
spot rates by maturity. If \\z_t\\ denotes the annual effective spot
rate for maturity \\t\\, then the discount factor is \\(1+z_t)^{-t}\\.

Computes \$\$ {}\_nE_x = (1+z_n)^{-n}\cdot {}\_np_x. \$\$

Computes \$\$ A\_{x:\overline{n}\|}^1 = \sum\_{t=1}^{n}(1+z_t)^{-t}\cdot
{}\_{t-1}p_x \cdot q\_{x+t-1}. \$\$

Computes \$\$ A\_{x:\overline{n}\|} = A\_{x:\overline{n}\|}^1 + {}\_nE_x
\$\$ using spot-rate discount factors.

For an immediate annuity, \$\$ a\_{x:\overline{n}\|} =
\sum\_{t=1}^{n}(1+z_t)^{-t}\cdot {}\_tp_x. \$\$

## Usage

``` r
nEx_spot(qx, z, benefit = 1)

Axn1_spot(qx, z, benefit = 1)

Axn_spot(qx, z, benefit = 1)

axn_spot(qx, z, type = c("immediate", "due"), benefit = 1)
```

## Arguments

- qx:

  Numeric vector of one-year mortality rates.

- z:

  Numeric vector of annual effective spot rates for maturities
  \\1,\dots,n\\.

- benefit:

  Amount of each annuity payment.

- type:

  Either `"immediate"` or `"due"`.

## Value

A numeric scalar.

A numeric scalar.

A numeric scalar.

A numeric scalar.

## Details

For an annuity-due, \$\$ \ddot{a}\_{x:\overline{n}\|} =
\sum\_{t=0}^{n-1}(1+z_t)^{-t}\cdot {}\_tp_x, \$\$ where the time-0
discount factor is 1.

## Examples

``` r
qx <- c(.02, .03, .04, .05, .06)
z  <- c(.03, .04, .05, .06, .07)
nEx_spot(qx, z, benefit = 1000000)
#> [1] 581034.1

qx <- c(.02, .03, .04, .05, .06)
z  <- c(.03, .04, .05, .06, .07)
Axn1_spot(qx, z)
#> [1] 0.1526756

qx <- c(.02, .03, .04, .05, .06)
z  <- c(.03, .04, .05, .06, .07)
Axn_spot(qx, z)
#> [1] 0.7337096

qx <- c(.02, .03, .04, .05, .06)
z  <- c(.03, .04, .05, .06, .07)
axn_spot(qx, z, type = "due")
#> [1] 4.30536
```

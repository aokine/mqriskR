# Variable-interest actuarial present value functions

Chapter 15 functions for actuarial present values under variable annual
effective interest rates interpreted as a yearly scenario \\i_1, i_2,
\dots, i_n\\.

Computes the APV of an \\n\\-year pure endowment under a variable annual
interest scenario: \$\$ {}\_nE_x = v_n \cdot {}\_np_x. \$\$

Computes the APV of an \\n\\-year term insurance with benefit paid at
the end of the year of death under a variable annual interest scenario:
\$\$ A\_{x:\overline{n}\|}^1 = \sum\_{t=1}^{n} v_t \cdot {}\_{t-1}p_x
\cdot q\_{x+t-1}. \$\$

Computes the APV of an \\n\\-year endowment insurance under a variable
annual interest scenario: \$\$ A\_{x:\overline{n}\|} =
A\_{x:\overline{n}\|}^1 + {}\_nE_x. \$\$

Computes the APV of an \\n\\-year temporary life annuity under a
variable annual interest scenario.

## Usage

``` r
nEx_var(qx, i, benefit = 1)

Axn1_var(qx, i, benefit = 1)

Axn_var(qx, i, benefit = 1)

axn_var(qx, i, type = c("immediate", "due"), benefit = 1)
```

## Arguments

- qx:

  Numeric vector of one-year mortality rates \\q_x, q\_{x+1}, \dots,
  q\_{x+n-1}\\.

- i:

  Numeric vector of annual effective interest rates \\i_1, i_2, \dots,
  i_n\\.

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

If a benefit amount is supplied, the function returns that benefit times
the APV factor.

If a benefit amount is supplied, the function returns that benefit times
the APV factor.

If a benefit amount is supplied, the function returns that benefit times
the APV factor.

For an immediate annuity, \$\$ a\_{x:\overline{n}\|} = \sum\_{t=1}^{n}
v_t \cdot {}\_tp_x. \$\$

For an annuity-due, \$\$ \ddot{a}\_{x:\overline{n}\|} =
\sum\_{t=0}^{n-1} v_t \cdot {}\_tp_x, \$\$ with \\v_0 = 1\\.

## Examples

``` r
qx <- c(.03, .04, .05, .06, .07)
nEx_var(qx = qx, i = c(.06, .07, .08, .09, .10), benefit = 1000)
#> [1] 526.5563

qx <- c(.03, .04, .05, .06, .07)
Axn1_var(qx = qx, i = c(.06, .07, .08, .09, .10))
#> [1] 0.1799082

qx <- c(.03, .04, .05, .06, .07)
Axn_var(qx = qx, i = c(.06, .07, .08, .09, .10))
#> [1] 0.7064644

qx <- rep(.02, 5)
axn_var(qx = qx, i = c(.06, .05, .04, .03, .03), type = "immediate")
#> [1] 4.110256
axn_var(qx = qx, i = c(.03, .04, .05, .06, .07), type = "due")
#> [1] 4.458454
```

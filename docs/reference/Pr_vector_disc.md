# Profit vector for a discrete profit-analysis model

Computes the Chapter 17 profit vector \$\$ \mathbf{Pr} = (Pr_0, Pr_1,
\dots, Pr_n) \$\$ where \\Pr_0\\ is the negative pre-contract expense
and the yearly expected profit values are calculated from the general
discrete expression in Equation (17.1).

## Usage

``` r
Pr_vector_disc(
  V,
  G,
  i,
  r = 0,
  e = 0,
  q1,
  q2 = 0,
  b1,
  b2 = 0,
  s1 = 0,
  s2 = 0,
  p_tau = NULL,
  pre_contract_expense = 0
)
```

## Arguments

- V:

  Vector of gross premium reserves \\{}\_tV^G\\ of length \\n+1\\,
  including the issue-time reserve and the terminal reserve.

- G:

  Gross premium vector for policy years 1 through \\n\\.

- i:

  Interest-rate vector for policy years 1 through \\n\\.

- r:

  Percent-of-premium expense vector.

- e:

  Fixed expense vector.

- q1:

  First decrement probabilities, typically death.

- q2:

  Second decrement probabilities, typically surrender or lapse. Defaults
  to 0.

- b1:

  Benefit vector for decrement 1.

- b2:

  Benefit vector for decrement 2. Defaults to 0.

- s1:

  Settlement-expense vector for decrement 1. Defaults to 0.

- s2:

  Settlement-expense vector for decrement 2. Defaults to 0.

- p_tau:

  Optional vector of in-force probabilities \\p\_{x+t}^{(\tau)}\\. If
  omitted, it is computed as \\1-q^{(1)}-q^{(2)}\\.

- pre_contract_expense:

  Positive pre-contract expense amount. The returned first element is
  \\Pr_0 = -\text{pre\\contract\\expense}\\.

## Value

Numeric vector of length \\n+1\\.

## Details

This implementation allows for two decrements, typically death and
withdrawal/surrender.

## Examples

``` r
V <- c(0, 5.66, 6.17, 0)
qx <- c(0.00142, 0.00153, 0.00166)
Pr_vector_disc(
  V = V, G = 95, i = 0.06, r = 0.05, e = 10,
  q1 = qx, b1 = 50000, pre_contract_expense = 15
)
#>        Pr0        Pr1        Pr2        Pr3 
#> -15.000000   8.413037   8.404040   8.605200 
```

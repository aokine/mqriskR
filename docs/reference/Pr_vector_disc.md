# Profit vector for a discrete profit-analysis model

Computes expected profit by policy year for a discrete contract with up
to two decrements. The first element is the negative pre-contract
expense.

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

  Numeric vector of gross premium reserves with length `n + 1`,
  including the issue-time and terminal reserves.

- G:

  Gross premium by policy year.

- i:

  Annual effective interest rate by policy year. Values must be greater
  than `-1`.

- r:

  Percent-of-premium expense rate by policy year. Values must lie in
  `[0, 1]`.

- e:

  Fixed expense by policy year.

- q1:

  Probability of the first decrement by policy year.

- q2:

  Probability of the second decrement by policy year.

- b1:

  Benefit payable on the first decrement.

- b2:

  Benefit payable on the second decrement.

- s1:

  Settlement expense associated with the first decrement.

- s2:

  Settlement expense associated with the second decrement.

- p_tau:

  Optional in-force probability by policy year. If omitted, it is
  calculated as `1 - q1 - q2`.

- pre_contract_expense:

  Nonnegative scalar pre-contract expense.

## Value

A named numeric vector of length `n + 1`.

## Details

For policy year \\k\\, the expected profit is \$\$ \[V\_{k-1} +
G_k(1-r_k)-e_k\](1+i_k) - \[(b_k^{(1)}+s_k^{(1)})q_k^{(1)}
+(b_k^{(2)}+s_k^{(2)})q_k^{(2)} +V_kp_k^{(\tau)}\]. \$\$

Scalar yearly inputs are recycled to the number of policy years
determined by `length(V) - 1`.

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

# Projected asset share path \\{}\_{k}AS\\

Computes projected asset shares recursively using Equations (14.5b) and
(14.6b) of Chapter 14, with optional support for a survival benefit
payable at the end of year \\k\\.

## Usage

``` r
AS_path(AS0, G, r, e, b1, b2, q1, q2, p_tau, i, b3 = NULL)
```

## Arguments

- AS0:

  Initial asset share \\{}\_{0}AS\\.

- G:

  Level annual premium.

- r:

  Numeric vector of percent-of-premium expense factors.

- e:

  Numeric vector of fixed contract expenses.

- b1:

  Numeric vector of Cause 1 benefit amounts.

- b2:

  Numeric vector of Cause 2 benefit amounts.

- q1:

  Numeric vector of Cause 1 decrement probabilities.

- q2:

  Numeric vector of Cause 2 decrement probabilities.

- p_tau:

  Numeric vector of in-force probabilities.

- i:

  Effective annual interest rate.

- b3:

  Optional numeric vector of survival benefit amounts payable at the end
  of year \\k\\ conditional on survival through year \\k\\. Defaults to
  a zero vector.

## Value

A data frame with columns `k` and `AS`.

## Details

For policy year \\k\\, \$\$ \[{}\_{k-1}AS + G(1-r_k) - e_k\](1+i) =
b_k^{(1)} q\_{x+k-1}^{(1)} + b_k^{(2)} q\_{x+k-1}^{(2)} +
p\_{x+k-1}^{(\tau)} \left(b_k^{(3)} + {}\_{k}AS\right) \$\$ so that \$\$
{}\_{k}AS = \frac{\[{}\_{k-1}AS + G(1-r_k) - e_k\](1+i) - b_k^{(1)}
q\_{x+k-1}^{(1)} - b_k^{(2)} q\_{x+k-1}^{(2)}}{p\_{x+k-1}^{(\tau)}} -
b_k^{(3)} \$\$

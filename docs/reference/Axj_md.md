# Discrete multiple-decrement insurance APV \\A\_{x}^{(j)}\\

Computes the actuarial present value of a benefit payable at the end of
the year of decrement if decrement occurs by Cause \\j\\, matching
Equation (14.3b) in Chapter 14.

## Usage

``` r
Axj_md(qj, ptau, i, benefit = 1)
```

## Arguments

- qj:

  Numeric vector of conditional probabilities \\q\_{x+k}^{(j)}\\ for
  Cause \\j\\.

- ptau:

  Numeric vector of survival probabilities \\{}\_{k}p\_{x}^{(\tau)}\\ of
  remaining in force to duration \\k\\.

- i:

  Effective annual interest rate.

- benefit:

  Benefit amount payable on decrement by Cause \\j\\.

## Value

A numeric scalar.

## Details

The function evaluates \$\$ A\_{x}^{(j)} = \sum\_{k=0}^{n-1} v^{k+1}
{}\_{k}p\_{x}^{(\tau)} q\_{x+k}^{(j)} \$\$ with an optional benefit
amount multiplier.

## Examples

``` r
q1 <- c(.02, .02, .02, .02, .02)
q2 <- c(.03, .04, .05, .06, .00)
q3 <- c(.00, .00, .00, .00, .98)
qtau <- q1 + q2 + q3

ptau <- numeric(length(qtau))
ptau[1] <- 1
for (k in 2:length(qtau)) {
  ptau[k] <- prod(1 - qtau[1:(k - 1)])
}

Axj_md(qj = q1, ptau = ptau, i = 0.06, benefit = 1000)
#> [1] 75.34884
```

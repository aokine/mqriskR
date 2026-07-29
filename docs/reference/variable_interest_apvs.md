# Actuarial present values under variable annual interest rates

Computes life-contingent actuarial present values using a specified
sequence of annual effective interest rates.

## Usage

``` r
nEx_var(qx, i, benefit = 1)

Axn1_var(qx, i, benefit = 1)

Axn_var(qx, i, benefit = 1)

axn_var(qx, i, type = c("immediate", "due"), benefit = 1)
```

## Arguments

- qx:

  Numeric vector of one-year mortality probabilities.

- i:

  Numeric vector of annual effective interest rates. Each value must be
  greater than `-1`.

- benefit:

  Nonnegative scalar benefit or annuity payment amount.

- type:

  Character string equal to `"immediate"` or `"due"`.

## Value

A numeric scalar.

## Details

The vectors `qx` and `i` represent one valuation scenario and must have
the same positive length.

`nEx_var()` computes a pure endowment.

`Axn1_var()` computes term insurance payable at the end of the year of
death.

`Axn_var()` computes endowment insurance.

`axn_var()` computes a temporary annuity-immediate or annuity-due.

Each year's payment is discounted using the cumulative product of the
annual effective interest rates supplied in `i`.

The pure endowment is \$\$ {}\_nE = {}\_np_x\\v_n, \$\$ where \\v_n\\ is
the cumulative discount factor implied by the sequence of annual
effective interest rates.

Term insurance is obtained by discounting each possible death benefit
using the cumulative discount factor applicable to its payment year.

Endowment insurance equals the sum of the corresponding term insurance
and pure endowment.

Temporary annuities discount each payment using the cumulative discount
factors derived from the interest-rate sequence.

## Examples

``` r
qx <- c(0.03, 0.04, 0.05, 0.06, 0.07)
rates <- c(0.06, 0.07, 0.08, 0.09, 0.10)

nEx_var(qx, rates, benefit = 1000)
#> [1] 526.5563
Axn1_var(qx, rates)
#> [1] 0.1799082
Axn_var(qx, rates)
#> [1] 0.7064644
axn_var(qx, rates, type = "due")
#> [1] 4.081115
```

# Actuarial present values under spot rates

Computes life-contingent actuarial present values using annual effective
spot rates by maturity.

## Usage

``` r
nEx_spot(qx, z, benefit = 1)

Axn1_spot(qx, z, benefit = 1)

Axn_spot(qx, z, benefit = 1)

axn_spot(qx, z, type = c("immediate", "due"), benefit = 1)
```

## Arguments

- qx:

  Numeric vector of one-year mortality probabilities.

- z:

  Numeric vector of annual effective spot rates for maturities
  `1, ..., n`. Each value must be greater than `-1`.

- benefit:

  Nonnegative scalar benefit or annuity payment amount.

- type:

  Character string equal to `"immediate"` or `"due"`.

## Value

A numeric scalar.

## Details

If \\z_t\\ denotes the annual effective spot rate for maturity \\t\\,
the corresponding discount factor is \$\$(1+z_t)^{-t}.\$\$

`nEx_spot()` computes a pure endowment.

`Axn1_spot()` computes term insurance payable at the end of the year of
death.

`Axn_spot()` computes endowment insurance.

`axn_spot()` computes a temporary annuity-immediate or annuity-due.

Each payment is discounted using the spot rate corresponding to its
maturity rather than a single level interest rate.

The pure endowment is \$\$ {}\_nE = {}\_np_x(1+z_n)^{-n}. \$\$

Term insurance is obtained by discounting each death benefit using the
spot rate corresponding to its payment year.

Endowment insurance equals the sum of the corresponding term insurance
and pure endowment.

Temporary annuities discount each payment using the spot rate for its
payment time.

## Examples

``` r
qx <- c(0.02, 0.03, 0.04, 0.05, 0.06)
spot <- c(0.03, 0.04, 0.05, 0.06, 0.07)

nEx_spot(qx, spot, benefit = 1000)
#> [1] 581.0341
Axn1_spot(qx, spot)
#> [1] 0.1526756
Axn_spot(qx, spot)
#> [1] 0.7337096
axn_spot(qx, spot, type = "due")
#> [1] 4.30536
```

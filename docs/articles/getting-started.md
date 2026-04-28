# Getting Started with mqriskR

## Overview

The `mqriskR` package provides functions for actuarial risk modeling,
including survival models, insurance and annuity values, premium
calculations, reserves, multiple-decrement models, and mortality
improvement projections.

The functions are designed to align with standard actuarial notation and
support teaching, exam preparation, and reproducible actuarial analysis.

This vignette provides a brief introduction to the main functionality of
the package.

For detailed derivations and additional examples, see *Models for
Quantifying Risk*.

## Survival Models

We begin by computing survival probabilities under a simple model.

``` r

library(mqriskR)

# Probability of surviving 10 years from age 40
tpx(10, x = 40, model = "uniform", omega = 100)
#> [1] 0.8333333
```

This gives the probability that a life aged 40 survives 10 years under
the uniform distribution of deaths with limiting age 100.

## Insurance Functions

We next compute insurance actuarial present values.

``` r

# Whole life insurance
Ax(40, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.3154882

# 10-year term insurance
Axn1(40, n = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.1286956
```

The first value is the actuarial present value of a whole life insurance
paying 1 at the end of the year of death, while the second value is the
present value of a 10-year term insurance that pays 1 if death occurs
within the term.

## Life Annuities

We now compute annuity values.

``` r

# Whole life annuity-immediate
ax(40, i = 0.05, model = "uniform", omega = 100)
#> [1] 13.37475

# 10-year temporary annuity
axn(40, n = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 7.065505
```

The first value represents the present value of a whole life
annuity-immediate, while the second gives the present value of a 10-year
temporary annuity-immediate.

## Premium Calculation

Using the equivalence principle, we compute a level annual net premium.

``` r

Ax_term <- Axn1(40, n = 10, i = 0.05, model = "uniform", omega = 100)
adotx_term <- adotxn(40, n = 10, i = 0.05, model = "uniform", omega = 100)

premium <- Ax_term / adotx_term
premium
#> [1] 0.01703695
```

This is the level annual net premium calculated using the equivalence
principle, where expected present value of benefits equals expected
present value of premiums. This corresponds to the identity

``` math
P_{x:\overline{n}|}^{1} = A_{x:\overline{n}|}^{1} / \ddot{a}_{x:\overline{n}|}.
```

## Reserve Calculation

We now compute a simple prospective reserve at time $`t = 5`$ for the
same contract.

``` r

t <- 5

# Prospective reserve: V_t = A_{x+t:n-t} - P * ä_{x+t:n-t}
Ax_future <- Axn1(40 + t, n = 10 - t, i = 0.05, model = "uniform", omega = 100)
adotx_future <- adotxn(40 + t, n = 10 - t, i = 0.05, model = "uniform", omega = 100)

V_t <- Ax_future - premium * adotx_future
V_t
#> [1] 0.003947702
```

This is the prospective reserve at time $`t = 5`$, representing the
expected future loss at that time based on remaining benefits and
premiums. This illustrates the standard prospective reserve formula:

``` math
{}_{t}V_{x:\overline{n}|}^{1} = A_{x+t:\overline{n-t}|}^{1} - P_{x:\overline{n}|}^{1} \ddot{a}_{x+t:\overline{n-t}|}.
```

## Multiple-Decrement Models

The package supports multiple-decrement tables.

``` r

x <- 45:50

qmat <- cbind(
  q1 = c(.011, .012, .013, .014, .015, .016),
  q2 = rep(0.10, 6)
)

tbl <- md_table(x, qmat, radix = 1000)

tbl
#>    x    q1  q2  qtau  ptau      ltau        d1        d2      dtau
#> 1 45 0.011 0.1 0.111 0.889 1000.0000 11.000000 100.00000 111.00000
#> 2 46 0.012 0.1 0.112 0.888  889.0000 10.668000  88.90000  99.56800
#> 3 47 0.013 0.1 0.113 0.887  789.4320 10.262616  78.94320  89.20582
#> 4 48 0.014 0.1 0.114 0.886  700.2262  9.803167  70.02262  79.82578
#> 5 49 0.015 0.1 0.115 0.885  620.4004  9.306006  62.04004  71.34605
#> 6 50 0.016 0.1 0.116 0.884  549.0544  8.784870  54.90544  63.69030
```

This table summarizes the multiple-decrement model, including
cause-specific decrement probabilities, total decrement probabilities,
survival probabilities, and the number of lives remaining at each age.

We can compute survival probabilities:

``` r

npxtau_md(tbl, x = 46, n = 3)
#> [1] 0.6978632
```

This gives the probability that a life aged 46 remains in force for 3
years when all causes of decrement are considered.

## Mortality Improvement

The package also includes mortality improvement projections.

``` r

qx_proj(
  qx_base = 0.02,
  AAx = 0.01,
  base_year = 2020,
  proj_year = 2030
)
#> [1] 0.01808764
```

This gives the projected one-year death probability after allowing for
mortality improvement between the base year and the projection year.

## Summary

The `mqriskR` package provides a unified framework for actuarial
modeling, covering survival models, insurance and annuity functions,
premium calculations, reserves, multiple-decrement models, and mortality
improvement.

Users can combine these functions to build more complex actuarial models
in a transparent and reproducible way.

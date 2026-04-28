# UDD multiplier for m-thly insurance approximations

Under UDD, \\A_x^{(m)} = (i / i^{(m)}) A_x\\ and similarly for term and
deferred insurance.

## Usage

``` r
udd_mthly_multiplier(i, m)
```

## Arguments

- i:

  Numeric vector of effective annual interest rates.

- m:

  Positive integer payment frequency.

## Value

Numeric vector equal to \\i / i^{(m)}\\.

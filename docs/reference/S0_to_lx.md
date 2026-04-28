# Convert survival probabilities to life-table values

Converts Chapter 5 survival function values \\S_0(x)\\ into Chapter 6
life-table values \\l_x = l_0 S_0(x)\\ using a chosen radix.

## Usage

``` r
S0_to_lx(S0, radix = 1e+05)
```

## Arguments

- S0:

  Numeric vector of survival probabilities.

- radix:

  Positive radix \\l_0\\.

## Value

Numeric vector of \\l_x\\ values.

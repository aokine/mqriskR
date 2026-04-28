# Convert life-table values to survival probabilities

Converts Chapter 6 life-table values \\l_x\\ into Chapter 5 survival
probabilities \\S_0(x) = l_x / l_0\\.

## Usage

``` r
lx_to_S0(lx)
```

## Arguments

- lx:

  Numeric vector of life-table survivor values.

## Value

Numeric vector of \\S_0(x)\\ values.

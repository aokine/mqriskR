# Convert between compound-interest quantities

Provides consistent conversions between: - effective interest rate i -
effective discount rate d - force of interest delta

## Usage

``` r
interest_convert(i = NULL, d = NULL, delta = NULL, m = NULL)
```

## Arguments

- i:

  Effective interest rate.

- d:

  Effective discount rate.

- delta:

  Force of interest.

- m:

  Optional compounding frequency for the nominal rate convertible
  m-thly.

## Value

A list with elements i, d, delta and, if m is supplied, im (the nominal
rate convertible m-thly).

## Details

Exactly one of i, d, or delta must be provided.

## Examples

``` r
interest_convert(i = 0.05)
#> $i
#> [1] 0.05
#> 
#> $d
#> [1] 0.04761905
#> 
#> $delta
#> [1] 0.04879016
#> 
interest_convert(d = 0.04761905)
#> $i
#> [1] 0.05
#> 
#> $d
#> [1] 0.04761905
#> 
#> $delta
#> [1] 0.04879017
#> 
interest_convert(delta = log(1.05))
#> $i
#> [1] 0.05
#> 
#> $d
#> [1] 0.04761905
#> 
#> $delta
#> [1] 0.04879016
#> 
```

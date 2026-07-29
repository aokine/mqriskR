# Convert between compound-interest quantities

Provides consistent conversions between effective interest rate,
effective discount rate, force of interest, and optional nominal
interest rate convertible m-thly.

## Usage

``` r
interest_convert(i = NULL, d = NULL, delta = NULL, m = NULL)
```

## Arguments

- i:

  Effective interest rate. May be scalar or vector.

- d:

  Effective discount rate. May be scalar or vector.

- delta:

  Force of interest. May be scalar or vector.

- m:

  Optional positive integer compounding frequency for the nominal rate
  convertible m-thly.

## Value

A list with elements \`i\`, \`d\`, \`delta\`, and, if \`m\` is supplied,
\`im\` and \`m\`.

## Details

Exactly one of \`i\`, \`d\`, or \`delta\` must be provided.

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
interest_convert(i = c(0.03, 0.05, 0.07))
#> $i
#> [1] 0.03 0.05 0.07
#> 
#> $d
#> [1] 0.02912621 0.04761905 0.06542056
#> 
#> $delta
#> [1] 0.02955880 0.04879016 0.06765865
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

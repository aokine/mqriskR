# Cost of insurance for Type B universal life

Computes the one-period cost of insurance for a Type B universal life
policy using Equation (16.5a): \$\$ COI_t = \frac{B q\_{x+t-1}}{1+i^q}.
\$\$

## Usage

``` r
coi_ul_typeB(B, qx, iq)
```

## Arguments

- B:

  Face amount.

- qx:

  Mortality rate for the period.

- iq:

  Interest rate used in the cost-of-insurance calculation.

## Value

Numeric vector.

## Examples

``` r
coi_ul_typeB(B = 100000, qx = 0.00076, iq = 0.03)
#> [1] 73.78641
```

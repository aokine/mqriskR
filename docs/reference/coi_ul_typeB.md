# Cost of insurance for Type B universal life

Computes the one-period cost of insurance for a Type B universal life
contract: \$\$ \mathrm{COI}\_t = \frac{B_t q\_{x+t-1}}{1+i_t^q}. \$\$

## Usage

``` r
coi_ul_typeB(B, qx, iq)
```

## Arguments

- B:

  Face amount.

- qx:

  Mortality probability for the period.

- iq:

  Interest rate used in the cost-of-insurance calculation. Values must
  be greater than `-1`.

## Value

A numeric vector of cost-of-insurance charges.

## Details

The arguments may be scalars or vectors. Scalar arguments are recycled
to the common length.

## Examples

``` r
coi_ul_typeB(B = 100000, qx = 0.00076, iq = 0.03)
#> [1] 73.78641
coi_ul_typeB(
  B = 100000,
  qx = c(0.00076, 0.00081),
  iq = c(0.03, 0.035)
)
#> [1] 73.78641 78.26087
```

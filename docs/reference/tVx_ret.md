# Retrospective whole life reserve

Computes the retrospective net level premium reserve \$\${}\_tV_x =
P_x\ddot{s}\_{x:\overline{t}\|} -
\frac{A\_{x:\overline{t}\|}^{1}}{{}\_tE_x}.\$\$

## Usage

``` r
tVx_ret(x, t, i, model = NULL, ..., tbl = NULL)
```

## Arguments

- x:

  Issue age. May be scalar or vector.

- t:

  Nonnegative integer duration. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the actuarial functions.

- tbl:

  Optional life table object. Supply by name.

## Value

Numeric vector of retrospective reserve values.

## Details

The mortality basis may be supplied through either a life table or a
parametric survival model.

## Examples

``` r
tVx_ret(
  40,
  t = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.07250474
```

# Multi-life pure endowments

Computes pure endowment actuarial present values for joint-life and
last-survivor statuses for two independent lives.

## Usage

``` r
nExy(x, y, n, i, tbl = NULL, model = NULL, ...)

nExybar(x, y, n, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

- n:

  Nonnegative integer term. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the survival model.

## Value

A numeric vector of actuarial present values.

## Details

`nExy()` computes an \\n\\-year joint-life pure endowment, payable at
time \\n\\ if both lives survive.

`nExybar()` computes an \\n\\-year last-survivor pure endowment, payable
at time \\n\\ if at least one life survives.

The joint-life pure endowment is \$\$ {}\_nE\_{xy} = v^n\\{}\_np\_{xy}.
\$\$

The last-survivor pure endowment is \$\$ {}\_nE\_{\overline{xy}} =
v^n\\{}\_np\_{\overline{xy}}. \$\$

Under independence, \$\$ {}\_np\_{xy} = {}\_np_x\\{}\_np_y \$\$ and \$\$
{}\_np\_{\overline{xy}} = {}\_np_x + {}\_np_y - {}\_np_x\\{}\_np_y. \$\$

## Examples

``` r
nExy(
  x = 40,
  y = 50,
  n = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.4092755

nExybar(
  x = 40,
  y = 50,
  n = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.5934495
```

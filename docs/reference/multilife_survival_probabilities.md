# Multi-life survival and failure probabilities

Computes joint-life and last-survivor survival and failure probabilities
for two independent lives.

## Usage

``` r
tpxy(x, y, t, tbl = NULL, model = NULL, ...)

tqxy(x, y, t, tbl = NULL, model = NULL, ...)

tpxybar(x, y, t, tbl = NULL, model = NULL, ...)

tqxybar(x, y, t, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

- t:

  Nonnegative duration. May be scalar or vector.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the survival model.

## Value

A numeric vector of probabilities.

## Details

`tpxy()` computes the joint-life survival probability \$\$ {}\_tp\_{xy}
= {}\_tp_x\\{}\_tp_y. \$\$

`tqxy()` computes the joint-life failure probability \$\$ {}\_tq\_{xy} =
1-{}\_tp\_{xy}. \$\$

`tpxybar()` computes the last-survivor survival probability \$\$
{}\_tp\_{\overline{xy}} = {}\_tp_x + {}\_tp_y - {}\_tp_x\\{}\_tp_y. \$\$

`tqxybar()` computes the last-survivor failure probability \$\$
{}\_tq\_{\overline{xy}} = 1-{}\_tp\_{\overline{xy}}. \$\$

All calculations assume the future lifetimes are independent.

Joint-life probabilities require both lives to satisfy the survival
condition, whereas last-survivor probabilities require at least one life
to survive.

## Examples

``` r
tpxy(
  40, 50,
  t = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.6666667

tqxy(
  40, 50,
  t = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.3333333

tpxybar(
  40, 50,
  t = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.9666667

tqxybar(
  40, 50,
  t = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.03333333
```

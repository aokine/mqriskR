# Continuous multi-life insurance

Computes continuous joint-life, last-survivor, and contingent whole-life
insurance values for two independent lives.

## Usage

``` r
Abarxy(x, y, i, tbl = NULL, model = NULL, ...)

Abarxybar(x, y, i, tbl = NULL, model = NULL, ...)

Abarxy1(x, y, i, tbl = NULL, model = NULL, ...)

Abaryx1(x, y, i, tbl = NULL, model = NULL, ...)

Abarxy2(x, y, i, tbl = NULL, model = NULL, ...)

Abaryx2(x, y, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object retained for backward compatibility.
  Continuous calculations currently require `model`.

- model:

  Parametric survival model.

- ...:

  Additional parameters passed to the survival and hazard functions.

## Value

A numeric vector of actuarial present values.

## Details

`Abarxy()` computes joint-life insurance payable at the first death.

`Abarxybar()` computes last-survivor insurance payable at the second
death.

`Abarxy1()` and `Abaryx1()` compute contingent insurance payable when
the specified life dies first.

`Abarxy2()` and `Abaryx2()` compute contingent insurance payable when
the specified life dies second.

These functions require a parametric survival model.

Under independence, `Abarxy()` represents insurance payable at the first
death and `Abarxybar()` represents insurance payable at the second
death.

The contingent functions distinguish both the life whose death triggers
payment and whether that death occurs first or second.

## Examples

``` r
Abarxy(
  x = 40,
  y = 50,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.4899262

Abarxy1(
  x = 40,
  y = 50,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.2137821
```

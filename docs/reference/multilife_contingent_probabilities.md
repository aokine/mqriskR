# Contingent multi-life probabilities

Computes probabilities associated with the order of death of two
independent lives over a specified term.

## Usage

``` r
tqxy1(x, y, n, tbl = NULL, model = NULL, ...)

tqyx1(x, y, n, tbl = NULL, model = NULL, ...)

tqxy2(x, y, n, tbl = NULL, model = NULL, ...)

tqyx2(x, y, n, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of the first life. May be scalar or vector.

- y:

  Age of the second life. May be scalar or vector.

- n:

  Term in years. May be scalar or vector.

- tbl:

  Optional life table object retained for backward compatibility.
  Continuous calculations currently require `model`.

- model:

  Parametric survival model.

- ...:

  Additional parameters passed to the survival and hazard functions.

## Value

A numeric vector of contingent probabilities.

## Details

`tqxy1()` computes the probability that the first life dies before the
second life within \\n\\ years.

`tqyx1()` computes the probability that the second life dies before the
first life within \\n\\ years.

`tqxy2()` computes the probability that the first life dies after the
second life but within \\n\\ years.

`tqyx2()` computes the probability that the second life dies after the
first life but within \\n\\ years.

These functions require a parametric survival model because an annual
life table does not uniquely determine the within-year order of death
without an additional interpolation assumption.

All calculations assume the two future lifetimes are independent.

The probabilities are obtained by integrating the joint survival
function together with the appropriate force of mortality.

The four functions partition the probability that one of the two lives
dies within the specified term according to the order of death.

## Examples

``` r
tqxy1(
  40, 50,
  n = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.15

tqyx1(
  40, 50,
  n = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.1833333

tqxy2(
  40, 50,
  n = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.01666667

tqyx2(
  40, 50,
  n = 10,
  model = "uniform",
  omega = 100
)
#> [1] 0.01666667
```

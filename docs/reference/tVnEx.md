# Pure endowment net level premium reserve

Computes the Chapter 10 prospective reserve for an n-year pure
endowment.

## Usage

``` r
tVnEx(x, n, t, i, model, ...)
```

## Arguments

- x:

  Issue age.

- n:

  Term in years.

- t:

  Duration.

- i:

  Effective annual interest rate.

- model:

  Survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

## Examples

``` r
tVnEx(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.3265297
```

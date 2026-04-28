# Endowment insurance reserve by retrospective method

Computes the Chapter 10 retrospective reserve
\\{}\_tV\_{x:\overline{n}\|} = P\_{x:\overline{n}\|}
\ddot{s}\_{x:\overline{t}\|} - {}\_tk_x\\ for \\t \le n\\.

## Usage

``` r
tVxn_ret(x, n, t, i, model, ...)
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
tVxn_ret(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.3448973
```

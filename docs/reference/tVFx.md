# Full preliminary term reserve for whole life insurance

Computes the FPT reserve for a whole life insurance. For whole life
insurance, \\{}\_1V^F = 0\\ and for \\t \ge 1\\, \\{}\_tV^F =
{}\_{t-1}V\_{x+1}^{NLP}\\.

Computes the FPT reserve for a whole life insurance. For whole life
insurance, \\{}\_1V^F = 0\\ and for \\t \ge 1\\, \\{}\_tV^F =
{}\_{t-1}V\_{x+1}^{NLP}\\.

## Usage

``` r
tVFx(x, t, i, tbl = NULL, model = NULL, ...)

tVFx(x, t, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age.

- t:

  Duration.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

Numeric vector.

## Examples

``` r
tVFx(40, t = 5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.02773584
tVFx(40, t = 5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.02773584
```

# Distribution functions for age-at-failure T0

Convenience functions for the CDF and density of \\T_0\\.

## Usage

``` r
F0(t, model, ...)

f0(t, model, ...)
```

## Arguments

- t:

  Numeric vector of times (\\t \ge 0\\).

- model:

  One of `"uniform"`, `"exponential"`, `"gompertz"`, `"makeham"`,
  `"weibull"`.

- ...:

  Model parameters:

  - uniform: `omega`

  - exponential: `lambda`

  - gompertz: `B`, `c`

  - makeham: `A`, `B`, `c`

  - weibull: `shape`, `scale`

## Value

Numeric vector. For `F0`: CDF values in \\\[0,1\]\\. For `f0`: density
values (\\\ge 0\\).

## Details

- `F0(t)` computes \\F_0(t)=Pr(T_0 \le t)=1-S_0(t)\\

- `f0(t)` computes \\f_0(t)=\frac{d}{dt}F_0(t)\\

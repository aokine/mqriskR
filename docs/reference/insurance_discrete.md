# Discrete insurance models (Chapter 7)

Discrete contingent payment / insurance functions matching Chapter 7
notation.

## Details

These functions handle:

- whole life insurance: \\A_x\\,

- term insurance: \\A\_{x:\overline{n}\|}^{1}\\,

- deferred insurance: \\{}\_{n\mid}A_x\\,

- pure endowment: \\{}\_nE_x\\,

- endowment insurance: \\A\_{x:\overline{n}\|}\\,

- second moments and variances.

The functions may be evaluated from either:

- a life table object via `tbl = ...`, or

- a parametric survival model via `model = ...` and additional
  parameters.

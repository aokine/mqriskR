# m-thly insurance models (Chapter 7)

Exact m-thly contingent payment / insurance functions matching Chapter 7
notation.

## Details

These functions handle:

- m-thly whole life insurance: \\A_x^{(m)}\\,

- m-thly term insurance: \\A\_{x:\overline{n}\|}^{1(m)}\\,

- m-thly deferred insurance: \\{}\_{n\mid}A_x^{(m)}\\,

- m-thly endowment insurance: \\A\_{x:\overline{n}\|}^{(m)}\\,

- second moments and variances.

These functions are evaluated exactly from the parametric survival-model
framework. For UDD approximations from annual life tables, use the
corresponding \*\_udd helpers in insurance_utils.R.

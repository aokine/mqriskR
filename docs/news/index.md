# Changelog

## mqriskR 0.1.1

CRAN release: 2026-07-27

### Improvements

#### Life Table Support

- Improved support for life table objects across many actuarial
  functions.
- Made behavior more consistent between life table and parametric
  survival model implementations.

#### Vectorization

- Added vectorized input support to many functions where appropriate.
- Standardized argument recycling across related functions while
  preserving existing scalar behavior.

#### Finite Life Tables

- Improved handling of boundary cases near the limiting age.
- Made insurance and annuity functions behave more consistently for
  finite life tables.

#### Documentation

- Reorganized documentation using descriptive actuarial topic names
  instead of chapter-based names.
- Updated shared documentation pages for improved readability.
- Improved mathematical notation, examples, and parameter descriptions.
- Removed obsolete documentation pages generated from previous
  documentation structures.

#### Code Quality

- Improved consistency across related actuarial functions.
- Reduced duplicated code where practical.
- Enhanced input validation and error handling.

#### Testing

- Expanded test coverage for edge cases and vectorized inputs.
- Improved consistency of package examples and documentation.

### Compatibility

- No breaking changes.
- Existing function names and syntax remain unchanged.
- Code written for previous versions of the package should continue to
  work.

# Changelog

All notable changes to **mqriskR** will be documented in this file.

The format is based on the principles of *Keep a Changelog* and versioning follows semantic versioning.

---

## [0.1.1] - 2026-07-27

### Improved

#### Life Table Support

- Improved support for life table objects across actuarial functions.
- Standardized behavior between life table and parametric survival model implementations.

#### Vectorization

- Added vectorized input support to functions where mathematically appropriate.
- Standardized argument recycling across related functions while preserving scalar behavior.

#### Finite Life Tables

- Improved handling of boundary cases near the limiting age.
- Made insurance, annuity, premium, and reserve calculations more consistent for finite life tables.

#### Documentation

- Replaced chapter-based documentation topics with descriptive actuarial topic names.
- Consolidated shared documentation pages for improved readability.
- Updated parameter descriptions, mathematical notation, and examples.
- Removed obsolete generated documentation pages.
- Improved the organization of the package manual.

#### Examples

- Simplified examples.
- Standardized parameter names across examples.
- Improved consistency between examples.

#### Code Quality

- Reduced duplicated code where practical.
- Improved consistency across related actuarial functions.
- Enhanced input validation and error handling.

#### Testing

- Expanded testing recommendations for vectorized inputs and edge cases.
- Improved consistency between tests and documented examples.

### Compatibility

- No breaking API changes.
- Existing exported function names remain unchanged.
- Existing syntax remains unchanged.
- Code written for previous versions continues to work.

---

## [0.1.0]

### Initial Release

- Initial public release of the **mqriskR** package.
- Core actuarial functions for life contingencies.
- Survival models.
- Life insurance and annuity calculations.
- Multiple-decrement models.
- Mortality improvement functions.
- Standard actuarial notation throughout.

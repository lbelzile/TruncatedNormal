## Version 2.3

### Bug fixes

- Remove ghost argument `log` (never implemented) from functions `pmvnorm` and `pmvt`; these now return a warning (thanks to Elysée Aristide)

## Version 2.2

### Bug fixes

- `ptmvnorm` and `ptmvt` now return 1 if all points exceed the upper bounds (incorrectly returned zero in v2.1.2 and `NA` in previous versions) (thanks to Adelchi Azzalini).
- Probability estimates now capped at 0/1 in ptmvt/ptmvnorm.
- Check for solution to saddlepoint problem; the code now considers convex constrained optimization program for normal samples (QMC or MC).

### Changes

- No more warning for univariate Student distribution function evaluation.
- Unit testing for wrappers.

## Version 2.1.2

### Bug fixes 

- `ptmvtnorm` and `ptmvt` return 0 for points outside the truncated region instead of `NA` (thanks to Adelchi Azzalini)
- Fixed a bug in `rtmvt` that returned incorrect output (tranposed mean vector)


## Version 2.1

### Bug fixes

- Cholesky decomposition with permutation now checks the arguments to ensure `lb` less than `ub`, to avoid segfault errors. In case of degenerate bounds, the conditional distribution is used with fixed components for the degenerate variables (issue #2)

## Changes 

- Moved back probit example in documentation
- Change back to `randtoolbox` package following its reuploading on **CRAN**

## Version 2.0
### Bug fixes

- Throw error if `NaN` in bounds
- Handle univariate cases in `mvrandt` and `mvrandn`

### Changes

- Can specify mean vector directly in the functions
- Functions for multivariate Student simulation ported from Matlab
- Change to the interface
- Many functions now internal (documented)
- Cholesky matrix permutation (GGB and GE reordering now supported)
- Selected code rewritten in Rcpp
- New vignette
- Removed random seed
- Change to qrng package for QRN generation of the Sobol sequence (`randtoolbox` package unmaintained), but the sequence is not scrambled
- Use `nleq` solver to solve convex program (also check KKT conditions).

## Version 1.0

- Initial release

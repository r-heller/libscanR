# Create a LIBS Calibration Model

Constructs a `libs_calibration` S3 object holding a fitted calibration
model plus the underlying training data and figures of merit.

## Usage

``` r
ls_calibration(
  element,
  wavelength_nm,
  concentrations,
  intensities,
  model,
  method = "univariate",
  lod = NULL,
  loq = NULL,
  r_squared = NULL,
  unit = "ppm"
)
```

## Arguments

- element:

  Character. Element symbol (e.g. "Ca", "Fe").

- wavelength_nm:

  Numeric. Emission line wavelength used.

- concentrations:

  Numeric vector. Known concentrations of standards.

- intensities:

  Numeric vector. Measured peak intensities or areas (or ratios for
  `internal_std`). Same length as `concentrations`.

- model:

  Fitted model object (e.g. `lm`, `mvr`).

- method:

  Character. One of `"univariate"`, `"internal_std"`, `"pls"`,
  `"cf_libs"`. Default `"univariate"`.

- lod:

  Numeric. Limit of detection (3-sigma). Default `NULL`.

- loq:

  Numeric. Limit of quantification (10-sigma). Default `NULL`.

- r_squared:

  Numeric. Coefficient of determination. Default `NULL`.

- unit:

  Character. Concentration unit, e.g. `"ppm"`, `"wt_pct"`, `"mg_kg"`.
  Default `"ppm"`.

## Value

A `libs_calibration` S3 object.

## Examples

``` r
conc <- c(100, 500, 1000, 2500, 5000)
ints <- conc * 0.12 + stats::rnorm(5, 0, 5)
m <- stats::lm(conc ~ ints)
cal <- ls_calibration("Ca", 393.37, conc, ints, model = m,
                      r_squared = summary(m)$r.squared)
print(cal)
#> <libs_calibration>
#> • Element: "Ca"
#> • Line: 393.37 nm
#> • Method: "univariate"
#> • Unit: "ppm"
#> • Standards: 5
#> • R-squared: 0.9995
#> • LOD: (not set)
#> • LOQ: (not set)
```

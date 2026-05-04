#' Create a LIBS Calibration Model
#'
#' Constructs a `libs_calibration` S3 object holding a fitted calibration
#' model plus the underlying training data and figures of merit.
#'
#' @param element Character. Element symbol (e.g. "Ca", "Fe").
#' @param wavelength_nm Numeric. Emission line wavelength used.
#' @param concentrations Numeric vector. Known concentrations of standards.
#' @param intensities Numeric vector. Measured peak intensities or areas
#'   (or ratios for `internal_std`). Same length as `concentrations`.
#' @param model Fitted model object (e.g. `lm`, `mvr`).
#' @param method Character. One of `"univariate"`, `"internal_std"`, `"pls"`,
#'   `"cf_libs"`. Default `"univariate"`.
#' @param lod Numeric. Limit of detection (3-sigma). Default `NULL`.
#' @param loq Numeric. Limit of quantification (10-sigma). Default `NULL`.
#' @param r_squared Numeric. Coefficient of determination. Default `NULL`.
#' @param unit Character. Concentration unit, e.g. `"ppm"`, `"wt_pct"`,
#'   `"mg_kg"`. Default `"ppm"`.
#'
#' @return A `libs_calibration` S3 object.
#'
#' @examples
#' conc <- c(100, 500, 1000, 2500, 5000)
#' ints <- conc * 0.12 + stats::rnorm(5, 0, 5)
#' m <- stats::lm(conc ~ ints)
#' cal <- ls_calibration("Ca", 393.37, conc, ints, model = m,
#'                       r_squared = summary(m)$r.squared)
#' print(cal)
#'
#' @export
ls_calibration <- function(element, wavelength_nm, concentrations, intensities,
                           model, method = "univariate", lod = NULL, loq = NULL,
                           r_squared = NULL, unit = "ppm") {
  method <- match.arg(method, c("univariate", "internal_std", "pls", "cf_libs"))
  if (!is.character(element) || length(element) != 1) {
    cli::cli_abort("{.arg element} must be a single character string.")
  }
  if (!is.numeric(wavelength_nm) || length(wavelength_nm) != 1) {
    cli::cli_abort("{.arg wavelength_nm} must be a single numeric value.")
  }
  if (length(concentrations) != length(intensities)) {
    cli::cli_abort("{.arg concentrations} and {.arg intensities} must have equal length.")
  }

  structure(
    list(
      element = element,
      wavelength_nm = wavelength_nm,
      concentrations = as.numeric(concentrations),
      intensities = as.numeric(intensities),
      model = model,
      method = method,
      lod = lod,
      loq = loq,
      r_squared = r_squared,
      unit = unit,
      n_standards = length(concentrations)
    ),
    class = "libs_calibration"
  )
}

#' Print a LIBS Calibration
#' @param x A `libs_calibration` object.
#' @param ... Unused.
#' @return Invisibly `x`.
#' @export
print.libs_calibration <- function(x, ...) {
  cli::cli_inform(c(
    "{.cls libs_calibration}",
    "*" = "Element: {.val {x$element}}",
    "*" = "Line: {.val {round(x$wavelength_nm, 3)}} nm",
    "*" = "Method: {.val {x$method}}",
    "*" = "Unit: {.val {x$unit}}",
    "*" = "Standards: {x$n_standards}",
    "*" = if (!is.null(x$r_squared)) "R-squared: {.val {round(x$r_squared, 4)}}" else "R-squared: (not set)",
    "*" = if (!is.null(x$lod)) "LOD: {.val {signif(x$lod, 3)}} {x$unit}" else "LOD: (not set)",
    "*" = if (!is.null(x$loq)) "LOQ: {.val {signif(x$loq, 3)}} {x$unit}" else "LOQ: (not set)"
  ))
  invisible(x)
}

#' Summary of a LIBS Calibration
#' @param object A `libs_calibration` object.
#' @param ... Unused.
#' @return Invisibly a list of key figures of merit.
#' @export
summary.libs_calibration <- function(object, ...) {
  print(object)
  if (!is.null(object$model) && inherits(object$model, "lm")) {
    coefs <- stats::coef(object$model)
    cli::cli_inform(c(
      "i" = "Model: concentration = {round(coefs[1], 4)} + {round(coefs[2], 6)} * intensity"
    ))
  }
  invisible(list(
    element = object$element,
    wavelength_nm = object$wavelength_nm,
    method = object$method,
    r_squared = object$r_squared,
    lod = object$lod,
    loq = object$loq,
    n_standards = object$n_standards
  ))
}

#' Predict Concentrations from a Calibration Model
#'
#' @param object A `libs_calibration` object.
#' @param newdata A numeric vector of new intensities (or a matrix for PLS).
#' @param ... Unused.
#' @return Numeric vector of predicted concentrations.
#' @export
predict.libs_calibration <- function(object, newdata, ...) {
  if (is.null(object$model)) {
    cli::cli_abort("Calibration has no fitted model.")
  }
  if (object$method %in% c("univariate", "internal_std")) {
    df <- data.frame(intensity = as.numeric(newdata))
    pred <- stats::predict(object$model, newdata = df)
    return(as.numeric(pred))
  }
  if (object$method == "pls") {
    .require_pkg("pls", "Needed to predict from a PLS calibration.")
    if (is.null(dim(newdata))) {
      newdata <- matrix(newdata, nrow = 1)
    }
    pred <- stats::predict(object$model, newdata = newdata,
                           ncomp = object$model$ncomp)
    return(as.numeric(pred))
  }
  cli::cli_abort("Unsupported calibration method: {.val {object$method}}")
}

#' Test whether an object is a libs_calibration
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_libs_calibration <- function(x) {
  inherits(x, "libs_calibration")
}

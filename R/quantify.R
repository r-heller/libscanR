# Quantification and figures of merit
# -----------------------------------------------------------------------------

#' Quantify Element Concentration
#'
#' Applies a calibration model to unknown spectra.
#'
#' @param calibration A [ls_calibration()] object.
#' @param x A [ls_spectrum()] or [ls_dataset()] object.
#' @param window_nm Numeric. Peak integration window. Defaults to 1 if NULL.
#'
#' @return A [tibble::tibble()] with columns `sample_id`, `element`,
#'   `concentration`, `unit`, `below_lod`, `below_loq`.
#'
#' @examples
#' ds <- ls_example_data("calibration")
#' conc <- ds$sample_info$concentration
#' cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
#' ls_quantify(cal, ds)
#' @export
ls_quantify <- function(calibration, x, window_nm = NULL) {
  if (!inherits(calibration, "libs_calibration")) {
    cli::cli_abort("{.arg calibration} must be a {.cls libs_calibration} object.")
  }
  if (is.null(window_nm)) window_nm <- 1

  if (inherits(x, "libs_spectrum")) {
    inten <- ls_peak_area(x, calibration$wavelength_nm, window_nm)
    pred <- predict(calibration, inten)
    sid <- as.character(.meta_get(x$metadata, "sample_id", "unknown"))
    return(tibble::tibble(
      sample_id = sid,
      element = calibration$element,
      concentration = pred,
      unit = calibration$unit,
      below_lod = !is.null(calibration$lod) && pred < calibration$lod,
      below_loq = !is.null(calibration$loq) && pred < calibration$loq
    ))
  }
  if (inherits(x, "libs_dataset")) {
    intens <- vapply(x$spectra, function(s) {
      ls_peak_area(s, calibration$wavelength_nm, window_nm)
    }, numeric(1))
    preds <- predict(calibration, intens)
    sids <- vapply(x$spectra, function(s) {
      as.character(.meta_get(s$metadata, "sample_id", NA_character_))
    }, character(1))
    return(tibble::tibble(
      sample_id = sids,
      element = calibration$element,
      concentration = preds,
      unit = calibration$unit,
      below_lod = !is.null(calibration$lod) & preds < calibration$lod,
      below_loq = !is.null(calibration$loq) & preds < calibration$loq
    ))
  }
  cli::cli_abort("{.arg x} must be a {.cls libs_spectrum} or {.cls libs_dataset}.")
}

#' Limit of Detection (3-sigma)
#'
#' Calculates LOD using the 3-sigma criterion.
#'
#' @param calibration A [ls_calibration()] object.
#' @param blank Optional [ls_spectrum()] of a blank sample. If `NULL` (the
#'   default), LOD is recomputed from the residual standard error of the
#'   calibration model.
#' @param window_nm Numeric. Integration window for the blank. Default 1.
#'
#' @return Numeric LOD in the calibration's concentration unit.
#'
#' @examples
#' ds <- ls_example_data("calibration")
#' conc <- ds$sample_info$concentration
#' cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
#' ls_lod(cal)
#' @export
ls_lod <- function(calibration, blank = NULL, window_nm = 1) {
  if (!inherits(calibration, "libs_calibration")) {
    cli::cli_abort("{.arg calibration} must be a {.cls libs_calibration}.")
  }
  if (calibration$method == "pls") {
    return(calibration$lod)
  }
  slope <- tryCatch(stats::coef(calibration$model)[["intensity"]],
                    error = function(e) NA_real_)
  if (is.na(slope) || slope == 0) return(NA_real_)

  sigma <- if (is.null(blank)) {
    summary(calibration$model)$sigma
  } else {
    .validate_spectrum(blank)
    y <- .mean_intensity(blank)
    wl <- blank$wavelength
    lo <- calibration$wavelength_nm - window_nm / 2
    hi <- calibration$wavelength_nm + window_nm / 2
    idx <- which(wl >= lo & wl <= hi)
    stats::sd(y[idx], na.rm = TRUE)
  }
  3 * sigma / abs(slope)
}

#' Limit of Quantification (10-sigma)
#'
#' Calculates LOQ using the 10-sigma criterion.
#'
#' @param calibration A [ls_calibration()] object.
#' @param blank Optional [ls_spectrum()] of a blank sample.
#' @param window_nm Numeric. Integration window. Default 1.
#'
#' @return Numeric LOQ in the calibration's concentration unit.
#'
#' @examples
#' ds <- ls_example_data("calibration")
#' conc <- ds$sample_info$concentration
#' cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
#' ls_loq(cal)
#' @export
ls_loq <- function(calibration, blank = NULL, window_nm = 1) {
  lod <- ls_lod(calibration, blank, window_nm)
  lod * 10 / 3
}

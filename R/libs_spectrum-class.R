#' Create a LIBS Spectrum Object
#'
#' Constructs a `libs_spectrum` S3 object from a wavelength axis and one or
#' more intensity vectors (shots).
#'
#' @param wavelength Numeric vector of wavelengths in nanometers. Must be
#'   sortable and contain no NA.
#' @param intensity Numeric vector (single shot) or matrix (multiple shots)
#'   of emission intensities. If matrix: rows = shots, columns = wavelength
#'   channels, `ncol(intensity)` must equal `length(wavelength)`.
#' @param metadata Named list of metadata (e.g. `sample_id`, `material`,
#'   `gate_delay_us`, `integration_time_us`, `laser_energy_mj`, `atmosphere`,
#'   `date`). Default `list()`.
#' @param baseline Optional numeric vector of baseline values corresponding
#'   to each wavelength channel (populated by [ls_baseline()]). Default `NULL`.
#'
#' @return An object of class `libs_spectrum`: a named list with elements
#'   `wavelength`, `intensity` (matrix: shots x channels), `metadata`,
#'   `baseline`, `n_shots`, `n_channels`, and `range_nm`.
#'
#' @examples
#' wl <- seq(200, 900, length.out = 512)
#' int <- exp(-((wl - 393)^2) / 2) + stats::rnorm(512, 0, 0.05)
#' spec <- ls_spectrum(wl, int, metadata = list(sample_id = "demo"))
#' print(spec)
#'
#' @export
ls_spectrum <- function(wavelength, intensity, metadata = list(),
                        baseline = NULL) {
  if (!is.numeric(wavelength)) {
    cli::cli_abort("{.arg wavelength} must be numeric.")
  }
  if (anyNA(wavelength)) {
    cli::cli_abort("{.arg wavelength} must not contain NA values.")
  }
  if (length(wavelength) < 2) {
    cli::cli_abort("{.arg wavelength} must have at least 2 channels.")
  }
  if (!is.numeric(intensity)) {
    cli::cli_abort("{.arg intensity} must be numeric.")
  }

  if (is.null(dim(intensity))) {
    if (length(intensity) != length(wavelength)) {
      cli::cli_abort(c(
        "Length of {.arg intensity} must equal length of {.arg wavelength}.",
        "i" = "Got {length(intensity)} intensities and {length(wavelength)} wavelengths."
      ))
    }
    intensity <- matrix(intensity, nrow = 1)
  } else {
    if (ncol(intensity) != length(wavelength)) {
      cli::cli_abort(c(
        "{.code ncol(intensity)} must equal {.code length(wavelength)}.",
        "i" = "Got ncol = {ncol(intensity)} and length = {length(wavelength)}."
      ))
    }
  }

  if (!is.null(baseline)) {
    if (length(baseline) != length(wavelength)) {
      cli::cli_abort("{.arg baseline} must have same length as {.arg wavelength}.")
    }
  }

  if (!is.list(metadata)) {
    cli::cli_abort("{.arg metadata} must be a list.")
  }

  sorted <- .ensure_sorted(as.numeric(wavelength), intensity)

  structure(
    list(
      wavelength = sorted$wavelength,
      intensity = sorted$intensity,
      metadata = metadata,
      baseline = baseline,
      n_shots = nrow(sorted$intensity),
      n_channels = length(sorted$wavelength),
      range_nm = range(sorted$wavelength)
    ),
    class = "libs_spectrum"
  )
}

#' Test whether an object is a libs_spectrum
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_libs_spectrum <- function(x) {
  inherits(x, "libs_spectrum")
}

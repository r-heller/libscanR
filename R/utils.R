# Internal helper functions (not exported)

#' Validate a libs_spectrum object
#' @param x Object to validate
#' @param arg Argument name for error message
#' @keywords internal
#' @noRd
.validate_spectrum <- function(x, arg = "x") {
  if (!inherits(x, "libs_spectrum")) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a {.cls libs_spectrum} object.",
      "x" = "You supplied a {.cls {class(x)[1]}}."
    ))
  }
  invisible(TRUE)
}

#' Validate a libs_dataset object
#' @keywords internal
#' @noRd
.validate_dataset <- function(x, arg = "x") {
  if (!inherits(x, "libs_dataset")) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a {.cls libs_dataset} object.",
      "x" = "You supplied a {.cls {class(x)[1]}}."
    ))
  }
  invisible(TRUE)
}

#' Convert a wavelength (nm) to an index on a wavelength axis
#' @param wavelength Numeric vector (wavelength axis)
#' @param nm Numeric scalar wavelength value
#' @return Integer index of closest channel
#' @keywords internal
#' @noRd
.wavelength_to_idx <- function(wavelength, nm) {
  which.min(abs(wavelength - nm))
}

#' Get intensity vector for a spectrum (mean across shots if multi-shot)
#' @keywords internal
#' @noRd
.mean_intensity <- function(x) {
  if (!is.matrix(x$intensity)) return(as.numeric(x$intensity))
  if (nrow(x$intensity) == 1) return(as.numeric(x$intensity[1, ]))
  colMeans(x$intensity)
}

#' Build intensity matrix for a list of spectra (one row per spectrum)
#' @keywords internal
#' @noRd
.build_intensity_matrix <- function(spectra) {
  mat <- do.call(rbind, lapply(spectra, .mean_intensity))
  rownames(mat) <- vapply(spectra, function(s) {
    id <- s$metadata$sample_id
    if (is.null(id)) NA_character_ else as.character(id)
  }, character(1))
  mat
}

#' Trapezoidal integration
#' @keywords internal
#' @noRd
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2) return(0)
  sum((x[-1] - x[-n]) * (y[-1] + y[-n]) / 2)
}

#' Estimate noise as median absolute deviation from a smoothed baseline
#' @keywords internal
#' @noRd
.estimate_noise <- function(y) {
  if (length(y) < 5) return(stats::sd(y, na.rm = TRUE))
  diffs <- diff(y)
  stats::mad(diffs, na.rm = TRUE) / sqrt(2)
}

#' Path to packaged extdata file
#' @keywords internal
#' @noRd
.extdata_path <- function(file) {
  system.file("extdata", file, package = "libscanR")
}

#' Ensure wavelength axis is sorted ascending
#' @keywords internal
#' @noRd
.ensure_sorted <- function(wavelength, intensity) {
  if (is.unsorted(wavelength)) {
    ord <- order(wavelength)
    wavelength <- wavelength[ord]
    if (is.matrix(intensity)) {
      intensity <- intensity[, ord, drop = FALSE]
    } else {
      intensity <- intensity[ord]
    }
  }
  list(wavelength = wavelength, intensity = intensity)
}

#' Safe requireNamespace with informative error
#' @keywords internal
#' @noRd
.require_pkg <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- c("Package {.pkg {pkg}} is required.")
    if (!is.null(reason)) msg <- c(msg, i = reason)
    msg <- c(msg, i = "Install it with {.code install.packages(\"{pkg}\")}.")
    cli::cli_abort(msg)
  }
  invisible(TRUE)
}

#' Compose a human-readable range string
#' @keywords internal
#' @noRd
.range_str <- function(x, digits = 1) {
  r <- range(x, na.rm = TRUE)
  sprintf("%.*f-%.*f", digits, r[1], digits, r[2])
}

#' Safe name getter from a metadata list
#' @keywords internal
#' @noRd
.meta_get <- function(meta, key, default = NA) {
  if (is.null(meta)) return(default)
  v <- meta[[key]]
  if (is.null(v)) default else v
}

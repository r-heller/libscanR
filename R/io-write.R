#' Write a LIBS Spectrum to CSV
#'
#' Exports a `libs_spectrum` to a CSV file with a `wavelength_nm` column plus
#' one column per shot.
#'
#' @param x A `libs_spectrum` object.
#' @param path Character. Output file path.
#' @param include_metadata Logical. Prepend metadata as comment lines.
#'   Default `TRUE`.
#'
#' @return Invisibly returns `path`.
#'
#' @examples
#' spec <- ls_simulate_spectrum(seed = 1)
#' tmp <- tempfile(fileext = ".csv")
#' ls_write_csv(spec, tmp)
#' unlink(tmp)
#'
#' @export
ls_write_csv <- function(x, path, include_metadata = TRUE) {
  .validate_spectrum(x)
  inten <- x$intensity
  df <- data.frame(wavelength_nm = x$wavelength)
  if (nrow(inten) == 1) {
    df$intensity <- as.numeric(inten[1, ])
  } else {
    for (i in seq_len(nrow(inten))) {
      df[[paste0("shot_", i)]] <- as.numeric(inten[i, ])
    }
  }

  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  if (include_metadata && length(x$metadata) > 0) {
    for (key in names(x$metadata)) {
      v <- x$metadata[[key]]
      if (length(v) == 1 && (is.atomic(v))) {
        writeLines(paste0("# ", key, ": ", v), con)
      }
    }
  }
  utils::write.csv(df, con, row.names = FALSE)
  invisible(path)
}

#' Export a LIBS Dataset to a Directory
#'
#' Writes each spectrum in a dataset as a separate CSV file, plus a
#' `sample_info.csv` metadata table.
#'
#' @param dataset A `libs_dataset` object.
#' @param dir Character. Output directory (created if missing).
#' @param overwrite Logical. Overwrite existing files. Default `FALSE`.
#'
#' @return Invisibly returns `dir`.
#'
#' @examples
#' ds <- ls_example_data("tissue")[1:3]
#' tmp <- tempfile()
#' ls_export_spectra(ds, tmp)
#' list.files(tmp)
#' unlink(tmp, recursive = TRUE)
#'
#' @export
ls_export_spectra <- function(dataset, dir, overwrite = FALSE) {
  .validate_dataset(dataset)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  for (i in seq_along(dataset$spectra)) {
    sp <- dataset$spectra[[i]]
    fid <- as.character(.meta_get(sp$metadata, "sample_id",
                                  paste0("spectrum_", i)))
    path <- file.path(dir, paste0(fid, ".csv"))
    if (file.exists(path) && !overwrite) {
      cli::cli_warn("File {.path {path}} exists; skipping (set {.arg overwrite = TRUE}).")
      next
    }
    ls_write_csv(sp, path)
  }
  utils::write.csv(dataset$sample_info,
                   file.path(dir, "sample_info.csv"),
                   row.names = FALSE)
  invisible(dir)
}

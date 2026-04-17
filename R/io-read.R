#' Read a LIBS Spectrum from File
#'
#' Reads a single spectrum file. Supports CSV, TSV, and TXT formats where
#' columns represent wavelength and intensity.
#'
#' @param path Character. Path to spectrum file.
#' @param format Character. One of `"auto"` (detect from extension), `"csv"`,
#'   `"tsv"`, `"txt"`. Default `"auto"`.
#' @param wavelength_col Column index or name for wavelength values. Default 1.
#' @param intensity_col Column index/indices or name(s) for intensity values.
#'   Multiple columns are interpreted as multiple shots. Default 2.
#' @param skip Integer. Header lines to skip before the data table. Default 0.
#' @param metadata Named list of additional metadata. Default `list()`.
#' @param verbose Logical. Emit import messages. Default `TRUE`.
#'
#' @return A `libs_spectrum` object.
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' wl <- seq(200, 900, length.out = 100)
#' int <- exp(-((wl - 393)^2) / 2)
#' utils::write.csv(data.frame(wavelength = wl, intensity = int),
#'                  tmp, row.names = FALSE)
#' spec <- ls_read_spectrum(tmp, verbose = FALSE)
#' unlink(tmp)
#'
#' @export
ls_read_spectrum <- function(path, format = "auto", wavelength_col = 1,
                             intensity_col = 2, skip = 0,
                             metadata = list(), verbose = TRUE) {
  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }
  format <- match.arg(format, c("auto", "csv", "tsv", "txt"))
  if (format == "auto") {
    ext <- tolower(tools::file_ext(path))
    format <- switch(ext,
                     csv = "csv",
                     tsv = "tsv",
                     txt = "txt",
                     "csv")
  }
  delim <- switch(format,
                  csv = ",",
                  tsv = "\t",
                  txt = NULL)

  if (is.null(delim)) {
    df <- tryCatch(
      utils::read.table(path, header = TRUE, sep = "", skip = skip,
                        stringsAsFactors = FALSE),
      error = function(e) utils::read.table(path, header = FALSE, sep = "",
                                            skip = skip,
                                            stringsAsFactors = FALSE)
    )
  } else {
    df <- tryCatch(
      readr::read_delim(path, delim = delim, skip = skip,
                        show_col_types = FALSE,
                        progress = FALSE),
      error = function(e) {
        utils::read.table(path, header = TRUE, sep = delim, skip = skip,
                          stringsAsFactors = FALSE)
      }
    )
  }
  df <- as.data.frame(df)

  wl <- .pick_col(df, wavelength_col)
  int_cols_df <- .pick_cols(df, intensity_col)
  int_cols_df <- as.data.frame(lapply(int_cols_df, as.numeric))
  if (ncol(int_cols_df) == 1) {
    intensity <- as.numeric(int_cols_df[[1]])
  } else {
    intensity <- t(as.matrix(int_cols_df))
    dimnames(intensity) <- NULL
  }

  meta <- utils::modifyList(list(
    source_file = path,
    import_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    sample_id = tools::file_path_sans_ext(basename(path))
  ), metadata)

  spec <- ls_spectrum(
    wavelength = as.numeric(wl),
    intensity = intensity,
    metadata = meta
  )

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Read {.path {basename(path)}}: {spec$n_channels} channels, {spec$n_shots} shot(s)."
    ))
  }
  spec
}

#' Read All Spectra from a Directory
#'
#' Batch-imports spectrum files from a directory into a `libs_dataset`.
#'
#' @param dir Character. Directory path.
#' @param pattern Character. File pattern (regex). Default `"\\.(csv|tsv|txt)$"`.
#' @param recursive Logical. Recurse into subdirectories. Default `FALSE`.
#' @param sample_info data.frame/tibble. Optional sample metadata table.
#'   Must have a `sample_id` column matching (or being a subset of) the
#'   filenames without extension. Default `NULL`.
#' @param ... Additional arguments passed to [ls_read_spectrum()].
#'
#' @return A `libs_dataset` object.
#'
#' @examples
#' dir <- tempfile()
#' dir.create(dir)
#' for (i in 1:3) {
#'   wl <- seq(200, 900, length.out = 100)
#'   int <- exp(-((wl - 393)^2) / 2) + stats::rnorm(100, 0, 0.01)
#'   utils::write.csv(data.frame(w = wl, i = int),
#'                    file.path(dir, paste0("s", i, ".csv")),
#'                    row.names = FALSE)
#' }
#' ds <- ls_read_dir(dir, verbose = FALSE)
#' unlink(dir, recursive = TRUE)
#'
#' @export
ls_read_dir <- function(dir, pattern = "\\.(csv|tsv|txt)$",
                        recursive = FALSE, sample_info = NULL, ...) {
  if (!dir.exists(dir)) {
    cli::cli_abort("Directory not found: {.path {dir}}")
  }
  files <- list.files(dir, pattern = pattern, full.names = TRUE,
                      recursive = recursive)
  if (length(files) == 0) {
    cli::cli_abort("No files matching {.val {pattern}} in {.path {dir}}")
  }

  spectra <- lapply(files, function(f) {
    ls_read_spectrum(f, ...)
  })

  ls_dataset(spectra, sample_info = sample_info)
}

#' Read a LIBS Spectrum from CSV (Convenience Wrapper)
#'
#' @param path Character. Path to CSV file.
#' @param ... Additional arguments passed to [ls_read_spectrum()].
#' @return A `libs_spectrum` object.
#' @export
ls_read_csv <- function(path, ...) {
  ls_read_spectrum(path, format = "csv", ...)
}

# Internal helpers

#' @keywords internal
#' @noRd
.pick_col <- function(df, col) {
  if (is.character(col)) {
    if (!col %in% names(df)) {
      cli::cli_abort("Column {.val {col}} not found.")
    }
    return(df[[col]])
  }
  if (is.numeric(col)) {
    if (col < 1 || col > ncol(df)) {
      cli::cli_abort("Column index {.val {col}} out of range.")
    }
    return(df[[col]])
  }
  cli::cli_abort("Column index must be character or numeric.")
}

#' @keywords internal
#' @noRd
.pick_cols <- function(df, cols) {
  if (is.character(cols)) {
    missing_cols <- setdiff(cols, names(df))
    if (length(missing_cols) > 0) {
      cli::cli_abort("Columns not found: {.val {missing_cols}}")
    }
    return(df[, cols, drop = FALSE])
  }
  if (is.numeric(cols)) {
    if (length(cols) == 1 && cols == 2 && ncol(df) > 2) {
      cols <- 2:ncol(df)
    }
    bad <- cols < 1 | cols > ncol(df)
    if (any(bad)) {
      cli::cli_abort("Column index out of range: {.val {cols[bad]}}")
    }
    return(df[, cols, drop = FALSE])
  }
  cli::cli_abort("Intensity column index must be character or numeric.")
}

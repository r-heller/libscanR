#' Launch the libscanR Shiny Application
#'
#' Starts an interactive Shiny application for exploring, preprocessing, and
#' analyzing LIBS spectral data. The app provides six tabs: Import,
#' Preprocessing, Peaks, Calibration, Chemometrics, and Export.
#'
#' @param data Optional. A [ls_spectrum()], [ls_dataset()], or file path to
#'   preload. Default `NULL` (start with upload interface).
#' @param port Integer. Port for the app. Default `NULL` (auto-select).
#' @param launch.browser Logical. Open browser on launch. Default `TRUE`.
#'
#' @return Invisible `NULL`. Launches a Shiny application as a side effect.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   ls_run_app()
#' }
#' }
#' @export
ls_run_app <- function(data = NULL, port = NULL, launch.browser = TRUE) {
  needed <- c("shiny", "bslib", "plotly", "DT")
  missing_pkgs <- needed[!vapply(needed, requireNamespace,
                                 logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    cli::cli_abort(c(
      "The Shiny app needs the following package(s) to run:",
      "i" = "{.pkg {missing_pkgs}}",
      "i" = "Install with: {.code install.packages(c({paste0('\"', missing_pkgs, '\"', collapse = ', ')}))}"
    ))
  }
  app_dir <- system.file("shiny", "libscanR", package = "libscanR")
  if (!nzchar(app_dir)) {
    cli::cli_abort(c(
      "Could not find the Shiny app directory.",
      "i" = "Try re-installing {.pkg libscanR}."
    ))
  }
  shiny::shinyOptions(libscanR.data = data)
  shiny::runApp(app_dir, port = port,
                launch.browser = launch.browser,
                display.mode = "normal")
  invisible(NULL)
}

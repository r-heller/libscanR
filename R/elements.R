#' NIST Emission Line Database
#'
#' Returns a curated subset of the NIST Atomic Spectra Database containing
#' analytically useful emission lines for LIBS analysis, with a focus on
#' biomedically relevant elements. Data compiled from the NIST Atomic Spectra
#' Database (<https://physics.nist.gov/asd>); values are for reference use
#' only and users requiring high-precision spectroscopy should consult NIST
#' ASD directly.
#'
#' @param elements Character vector. Restrict to specific element symbols
#'   (e.g. `c("Ca", "Fe")`). Default `NULL` (all).
#' @param range_nm Numeric vector of length 2. Wavelength range filter in nm.
#'   Default `c(190, 900)`.
#' @param min_aki Numeric. Minimum transition probability (10^8 s^-1).
#'   Default 0 (no filter).
#' @param ionization Integer vector of ionization states to include:
#'   1 = neutral, 2 = singly ionized. Default `c(1L, 2L)`.
#' @param persistent_only Logical. Keep only persistent lines. Default `FALSE`.
#'
#' @return A [tibble::tibble()] with columns `element`, `ionization`,
#'   `wavelength_nm`, `aki`, `ei_ev`, `ek_ev`, `persistent`.
#'
#' @examples
#' head(ls_element_db())
#' ls_element_db(elements = "Ca", persistent_only = TRUE)
#' ls_element_db(elements = c("Na", "K"), range_nm = c(580, 780))
#'
#' @source NIST Atomic Spectra Database: <https://physics.nist.gov/asd>
#' @export
ls_element_db <- function(elements = NULL, range_nm = c(190, 900),
                          min_aki = 0, ionization = c(1L, 2L),
                          persistent_only = FALSE) {
  path <- .extdata_path("nist_lines.rds")
  if (path == "" || !file.exists(path)) {
    cli::cli_abort("Internal NIST line database not found.")
  }
  db <- readRDS(path)
  if (!is.null(elements)) {
    db <- db[db$element %in% elements, ]
  }
  if (!is.null(range_nm)) {
    db <- db[db$wavelength_nm >= range_nm[1] & db$wavelength_nm <= range_nm[2], ]
  }
  if (min_aki > 0) {
    db <- db[db$aki >= min_aki, ]
  }
  if (!is.null(ionization)) {
    db <- db[db$ionization %in% as.integer(ionization), ]
  }
  if (persistent_only) {
    db <- db[isTRUE_vec(db$persistent), ]
  }
  tibble::as_tibble(db)
}

#' @keywords internal
#' @noRd
isTRUE_vec <- function(x) {
  vapply(x, isTRUE, logical(1))
}

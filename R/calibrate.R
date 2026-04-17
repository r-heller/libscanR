# Calibration curve construction (univariate, internal std, PLS, CF-LIBS)
# -----------------------------------------------------------------------------

#' Build Calibration Curve
#'
#' Constructs a calibration model from standards with known concentrations.
#'
#' @param dataset A [ls_dataset()] containing standard spectra. Spectra
#'   order must correspond to `concentrations`.
#' @param element Character. Target element symbol.
#' @param line_nm Numeric. Emission line wavelength (nm).
#' @param concentrations Numeric vector of known concentrations, one per
#'   spectrum in `dataset`.
#' @param method Character. `"univariate"` (default), `"internal_std"`,
#'   `"pls"`, or `"cf_libs"`.
#' @param internal_std_nm Numeric. Internal standard line wavelength,
#'   required for `method = "internal_std"`.
#' @param window_nm Numeric. Peak integration window. Default 1.
#' @param pls_window_nm Numeric vector of length 2. Spectral window for PLS.
#'   Default `c(line_nm - 5, line_nm + 5)`.
#' @param n_components Integer. Number of PLS components. Default 5.
#' @param unit Character. Concentration unit. Default `"ppm"`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A [ls_calibration()] object.
#'
#' @examples
#' ds <- ls_example_data("calibration")
#' conc <- ds$sample_info$concentration
#' cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
#' print(cal)
#' @export
ls_calibrate <- function(dataset, element, line_nm, concentrations,
                         method = "univariate", internal_std_nm = NULL,
                         window_nm = 1,
                         pls_window_nm = NULL, n_components = 5,
                         unit = "ppm", verbose = TRUE) {
  .validate_dataset(dataset)
  method <- match.arg(method, c("univariate", "internal_std", "pls", "cf_libs"))
  if (length(concentrations) != dataset$n_spectra) {
    cli::cli_abort("{.arg concentrations} length ({length(concentrations)}) must match number of spectra ({dataset$n_spectra}).")
  }

  if (method == "cf_libs") {
    cli::cli_abort("Use {.fn ls_saha_boltzmann} for CF-LIBS on a single spectrum.")
  }

  # Peak intensities
  intensities <- vapply(dataset$spectra, function(s) {
    ls_peak_area(s, line_nm, window_nm)
  }, numeric(1))

  if (method == "internal_std") {
    if (is.null(internal_std_nm)) {
      cli::cli_abort("{.arg internal_std_nm} is required for internal standard calibration.")
    }
    std_int <- vapply(dataset$spectra, function(s) {
      ls_peak_area(s, internal_std_nm, window_nm)
    }, numeric(1))
    std_int[std_int == 0] <- NA
    intensities <- intensities / std_int
  }

  if (method == "pls") {
    .require_pkg("pls", "Required for PLS calibration.")
    if (is.null(pls_window_nm)) {
      pls_window_nm <- c(line_nm - 5, line_nm + 5)
    }
    idx <- which(dataset$wavelength >= pls_window_nm[1] &
                   dataset$wavelength <= pls_window_nm[2])
    X <- dataset$intensity_matrix[, idx, drop = FALSE]
    n_components <- min(n_components, nrow(X) - 1, ncol(X))
    model_df <- data.frame(Y = as.numeric(concentrations))
    model_df$X <- X
    fit <- pls::plsr(Y ~ X, ncomp = n_components, data = model_df,
                     validation = "none")
    pred <- as.numeric(stats::predict(fit, ncomp = n_components))
    r2 <- .r_squared(concentrations, pred)
    lod <- stats::sd(stats::residuals(fit)[, , n_components], na.rm = TRUE) * 3
    loq <- lod * 10 / 3
    if (verbose) {
      cli::cli_inform(c("v" = "PLS calibration fit: ncomp = {n_components}, R2 = {round(r2, 4)}"))
    }
    return(ls_calibration(
      element = element,
      wavelength_nm = line_nm,
      concentrations = concentrations,
      intensities = rowMeans(X),
      model = fit,
      method = "pls",
      lod = lod,
      loq = loq,
      r_squared = r2,
      unit = unit
    ))
  }

  # Univariate / internal_std: linear regression of concentration on intensity
  df <- data.frame(concentration = as.numeric(concentrations),
                   intensity = as.numeric(intensities))
  df <- df[stats::complete.cases(df), ]
  fit <- stats::lm(concentration ~ intensity, data = df)
  r2 <- summary(fit)$r.squared
  slope <- stats::coef(fit)[["intensity"]]
  sigma <- summary(fit)$sigma
  lod <- 3 * sigma / abs(slope)
  loq <- 10 * sigma / abs(slope)

  if (verbose) {
    cli::cli_inform(c(
      "v" = "{method} calibration for {element}: R2 = {round(r2, 4)}, LOD = {signif(lod, 3)} {unit}, LOQ = {signif(loq, 3)} {unit}"
    ))
  }

  ls_calibration(
    element = element,
    wavelength_nm = line_nm,
    concentrations = concentrations,
    intensities = intensities,
    model = fit,
    method = method,
    lod = lod,
    loq = loq,
    r_squared = r2,
    unit = unit
  )
}

#' Calibration-Free LIBS via Saha-Boltzmann
#'
#' Implements the CF-LIBS method for semi-quantitative analysis without
#' reference standards. Uses the Boltzmann plot method to estimate plasma
#' temperature and then relative elemental concentrations from line
#' intensities. Electron density is provided as an input or assumed.
#'
#' @param x A [ls_spectrum()] object (ideally baseline-corrected).
#' @param elements Character vector. Elements to quantify.
#' @param lines_nm Named list. Emission lines per element for the Boltzmann
#'   plot, e.g. `list(Ca = c(393.37, 396.85, 422.67))`.
#' @param electron_density Numeric. Assumed electron density (cm^-3).
#'   Default `1e17`.
#' @param window_nm Numeric. Integration window per line. Default 0.5.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A [tibble::tibble()] with columns `element`, `temperature_k`,
#'   `concentration_rel` (relative mole fraction, sum normalized to 1),
#'   `electron_density`, `n_lines` (lines used), `r_squared` (Boltzmann fit).
#'
#' @examples
#' spec <- ls_simulate_spectrum(
#'   elements = c(Ca = 5000, Fe = 200, Na = 1000), seed = 2)
#' spec <- ls_baseline(spec)
#' ls_saha_boltzmann(spec,
#'   elements = c("Ca", "Fe"),
#'   lines_nm = list(Ca = c(422.673, 445.478, 487.813),
#'                   Fe = c(371.994, 404.581, 438.354)))
#' @export
ls_saha_boltzmann <- function(x, elements, lines_nm,
                              electron_density = 1e17,
                              window_nm = 0.5, verbose = TRUE) {
  .validate_spectrum(x)
  if (is.null(names(lines_nm)) || !all(elements %in% names(lines_nm))) {
    cli::cli_abort("{.arg lines_nm} must be a named list with an entry per element.")
  }

  db <- ls_element_db(elements = elements)

  results <- list()
  for (el in elements) {
    ls <- lines_nm[[el]]
    el_db <- db[db$element == el & db$ionization == 1L, ]
    if (nrow(el_db) == 0) {
      if (verbose) cli::cli_warn("No neutral lines in DB for {el}; skipping.")
      next
    }
    # For each requested line: find closest DB entry, compute ln(I*lambda/(g*A))
    pts <- list()
    for (wn in ls) {
      matched <- el_db[which.min(abs(el_db$wavelength_nm - wn)), ]
      if (nrow(matched) == 0 || abs(matched$wavelength_nm - wn) > 0.5) next
      intensity <- ls_peak_area(x, wn, window_nm)
      if (intensity <= 0) next
      lambda_m <- matched$wavelength_nm * 1e-9
      # g upper level assumed 2J+1 (defaulted to 1)
      g <- 1
      A <- matched$aki * 1e8  # s^-1
      E_upper <- matched$ek_ev  # eV
      y <- log(intensity * lambda_m / (g * A))
      pts[[length(pts) + 1]] <- data.frame(E = E_upper, y = y)
    }
    if (length(pts) < 2) {
      if (verbose) cli::cli_warn("Fewer than 2 usable lines for {el}; skipping.")
      next
    }
    bp <- do.call(rbind, pts)
    fit <- stats::lm(y ~ E, data = bp)
    slope <- stats::coef(fit)[["E"]]
    intercept <- stats::coef(fit)[["(Intercept)"]]
    # Slope = -1 / (k_B * T) where k_B in eV/K
    k_B_eV <- 8.617333262e-5
    T_k <- -1 / (slope * k_B_eV)
    r2 <- summary(fit)$r.squared
    results[[el]] <- list(
      element = el,
      temperature_k = T_k,
      concentration_raw = exp(intercept),
      n_lines = nrow(bp),
      r_squared = r2
    )
  }

  if (length(results) == 0) {
    cli::cli_abort("No element quantified; Boltzmann plot failed for all species.")
  }

  raw <- vapply(results, function(r) r$concentration_raw, numeric(1))
  rel <- raw / sum(raw)
  out <- tibble::tibble(
    element = vapply(results, function(r) r$element, character(1)),
    temperature_k = vapply(results, function(r) r$temperature_k, numeric(1)),
    concentration_rel = rel,
    electron_density = electron_density,
    n_lines = vapply(results, function(r) as.integer(r$n_lines), integer(1)),
    r_squared = vapply(results, function(r) r$r_squared, numeric(1))
  )
  out
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.r_squared <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

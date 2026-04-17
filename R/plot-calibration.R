# Calibration plots
# -----------------------------------------------------------------------------

#' Plot a Calibration Curve
#'
#' @param calibration A [ls_calibration()] object (univariate/internal_std).
#' @param show_lod Logical. Show LOD as a horizontal dashed line. Default `TRUE`.
#' @param show_loq Logical. Show LOQ as a horizontal dashed line. Default `TRUE`.
#' @param show_prediction_interval Logical. Show 95% prediction band.
#'   Default `TRUE`.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' ds <- ls_example_data("calibration")
#' conc <- ds$sample_info$concentration
#' cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
#' ls_plot_calibration(cal)
#' @export
ls_plot_calibration <- function(calibration, show_lod = TRUE,
                                show_loq = TRUE,
                                show_prediction_interval = TRUE) {
  if (!inherits(calibration, "libs_calibration")) {
    cli::cli_abort("{.arg calibration} must be a {.cls libs_calibration}.")
  }
  if (calibration$method == "pls") {
    cli::cli_warn("Calibration plot for PLS shows intensity means; use {.fn ls_plot_residuals} for diagnostics.")
  }

  df <- data.frame(
    intensity = calibration$intensities,
    concentration = calibration$concentrations
  )
  df <- df[stats::complete.cases(df), ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$intensity,
                                        y = .data$concentration)) +
    ggplot2::geom_point(size = 2.5, colour = "#2E86AB", alpha = 0.8) +
    ggplot2::labs(
      x = "Peak intensity (a.u.)",
      y = paste0("Concentration (", calibration$unit, ")"),
      title = paste0(calibration$element, " @ ",
                     round(calibration$wavelength_nm, 2), " nm"),
      subtitle = sprintf("R^2 = %.4f | LOD = %.3g | LOQ = %.3g %s",
                         calibration$r_squared %||% NA_real_,
                         calibration$lod %||% NA_real_,
                         calibration$loq %||% NA_real_,
                         calibration$unit)
    ) +
    theme_libs()

  if (calibration$method %in% c("univariate", "internal_std") &&
      !is.null(calibration$model)) {
    if (show_prediction_interval) {
      p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                                    colour = "#A23B72", fill = "#A23B72",
                                    alpha = 0.15, level = 0.95, se = TRUE)
    } else {
      p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                                    colour = "#A23B72", se = FALSE)
    }
  }

  if (show_lod && !is.null(calibration$lod)) {
    p <- p + ggplot2::geom_hline(yintercept = calibration$lod,
                                 linetype = "dashed", colour = "orange") +
      ggplot2::annotate("text",
                        x = min(df$intensity, na.rm = TRUE),
                        y = calibration$lod,
                        label = "LOD", vjust = -0.5, hjust = 0,
                        colour = "orange", size = 3)
  }
  if (show_loq && !is.null(calibration$loq)) {
    p <- p + ggplot2::geom_hline(yintercept = calibration$loq,
                                 linetype = "dashed", colour = "red") +
      ggplot2::annotate("text",
                        x = min(df$intensity, na.rm = TRUE),
                        y = calibration$loq,
                        label = "LOQ", vjust = -0.5, hjust = 0,
                        colour = "red", size = 3)
  }
  p
}

#' Plot Calibration Residuals
#'
#' Standard regression diagnostic: residuals vs fitted.
#'
#' @param calibration A [ls_calibration()] object.
#' @return A ggplot2 object.
#'
#' @examples
#' ds <- ls_example_data("calibration")
#' conc <- ds$sample_info$concentration
#' cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
#' ls_plot_residuals(cal)
#' @export
ls_plot_residuals <- function(calibration) {
  if (!inherits(calibration, "libs_calibration")) {
    cli::cli_abort("{.arg calibration} must be a {.cls libs_calibration}.")
  }
  if (calibration$method %in% c("univariate", "internal_std")) {
    fitted <- stats::fitted(calibration$model)
    resid <- stats::residuals(calibration$model)
  } else if (calibration$method == "pls") {
    ncomp <- calibration$model$ncomp
    fitted <- stats::fitted(calibration$model)[, , ncomp]
    resid <- calibration$concentrations - fitted
  } else {
    cli::cli_abort("Unsupported calibration method: {.val {calibration$method}}")
  }
  df <- data.frame(fitted_value = fitted, residual = resid)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted_value,
                                   y = .data$residual)) +
    ggplot2::geom_point(colour = "#2E86AB", size = 2, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey30") +
    ggplot2::labs(x = "Fitted", y = "Residual",
                  title = paste0("Residuals: ", calibration$element)) +
    theme_libs()
}

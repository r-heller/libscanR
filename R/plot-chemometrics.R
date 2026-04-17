# Chemometric plots (PCA, PLS-DA)
# -----------------------------------------------------------------------------

#' Plot PCA Scores
#'
#' @param pca A `libs_pca` object from [ls_pca()].
#' @param pc_x Integer. PC for x-axis. Default 1.
#' @param pc_y Integer. PC for y-axis. Default 2.
#' @param color_by Character. Grouping variable. Default `NULL`.
#' @param ellipses Logical. Draw 95% confidence ellipses per group. Default `TRUE`.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' pca <- ls_pca(ds, n_components = 4)
#' ls_plot_pca(pca, color_by = "tissue")
#' @export
ls_plot_pca <- function(pca, pc_x = 1, pc_y = 2, color_by = NULL,
                        ellipses = TRUE) {
  if (!inherits(pca, "libs_pca")) {
    cli::cli_abort("{.arg pca} must be a {.cls libs_pca} object.")
  }
  if (pc_x > pca$n_components || pc_y > pca$n_components) {
    cli::cli_abort("pc_x / pc_y exceed number of computed components ({pca$n_components}).")
  }
  scores <- pca$scores
  df <- data.frame(
    PC_x = scores[, pc_x],
    PC_y = scores[, pc_y]
  )
  var_x <- round(pca$variance_explained[pc_x] * 100, 1)
  var_y <- round(pca$variance_explained[pc_y] * 100, 1)

  if (!is.null(color_by)) {
    if (!color_by %in% names(pca$sample_info)) {
      cli::cli_abort("{.arg color_by} column {.val {color_by}} not found.")
    }
    df$group <- as.character(pca$sample_info[[color_by]])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC_x, y = .data$PC_y,
                                          colour = .data$group))
    if (ellipses) {
      p <- p + ggplot2::stat_ellipse(level = 0.95, alpha = 0.5)
    }
    p <- p + ggplot2::geom_point(size = 2.5, alpha = 0.8) +
      ggplot2::scale_colour_brewer(palette = "Set1", name = color_by)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC_x, y = .data$PC_y)) +
      ggplot2::geom_point(size = 2.5, colour = "#2E86AB", alpha = 0.8)
  }
  p +
    ggplot2::labs(
      x = sprintf("PC%d (%.1f%%)", pc_x, var_x),
      y = sprintf("PC%d (%.1f%%)", pc_y, var_y),
      title = "PCA Score Plot"
    ) +
    theme_libs()
}

#' Plot PCA Loadings
#'
#' Plots the loadings for a selected principal component across the full
#' wavelength axis, optionally annotating the top contributing channels.
#'
#' @param pca A `libs_pca` object.
#' @param pc Integer. Which PC's loadings to show. Default 1.
#' @param n_top Integer. Number of top-absolute contributors to annotate.
#'   Default 10.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' pca <- ls_pca(ds, n_components = 4)
#' ls_plot_loadings(pca, pc = 1, n_top = 5)
#' @export
ls_plot_loadings <- function(pca, pc = 1, n_top = 10) {
  if (!inherits(pca, "libs_pca")) {
    cli::cli_abort("{.arg pca} must be a {.cls libs_pca} object.")
  }
  loadings <- pca$loadings[, pc]
  df <- data.frame(
    wavelength_nm = pca$wavelength,
    loading = loadings
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$wavelength_nm,
                                        y = .data$loading)) +
    ggplot2::geom_line(colour = "#2E86AB", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey30") +
    ggplot2::labs(x = "Wavelength (nm)", y = paste0("PC", pc, " loading"),
                  title = sprintf("PC%d loadings (%.1f%% variance)",
                                  pc, pca$variance_explained[pc] * 100)) +
    theme_libs()

  if (n_top > 0) {
    top_idx <- order(abs(loadings), decreasing = TRUE)[seq_len(min(n_top,
                                                                   length(loadings)))]
    top_df <- df[top_idx, ]
    p <- p + ggplot2::geom_point(data = top_df,
                                 ggplot2::aes(x = .data$wavelength_nm,
                                              y = .data$loading),
                                 colour = "#A23B72", size = 2) +
      ggplot2::geom_text(data = top_df,
                         ggplot2::aes(x = .data$wavelength_nm,
                                      y = .data$loading,
                                      label = round(.data$wavelength_nm, 1)),
                         size = 3, vjust = -0.8, colour = "#A23B72")
  }
  p
}

#' Plot PCA Variance (Scree Plot)
#'
#' @param pca A `libs_pca` object.
#' @return A ggplot2 object.
#' @examples
#' ds <- ls_example_data("tissue")
#' pca <- ls_pca(ds, n_components = 6)
#' ls_plot_scree(pca)
#' @export
ls_plot_scree <- function(pca) {
  if (!inherits(pca, "libs_pca")) {
    cli::cli_abort("{.arg pca} must be a {.cls libs_pca} object.")
  }
  df <- data.frame(
    component = factor(paste0("PC", seq_along(pca$variance_explained)),
                       levels = paste0("PC", seq_along(pca$variance_explained))),
    variance_explained = pca$variance_explained * 100,
    cumulative = pca$cumulative_variance * 100
  )
  ggplot2::ggplot(df, ggplot2::aes(x = .data$component,
                                   y = .data$variance_explained)) +
    ggplot2::geom_col(fill = "#2E86AB", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = .data$cumulative, group = 1),
                       colour = "#A23B72", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = .data$cumulative),
                        colour = "#A23B72", size = 2) +
    ggplot2::labs(x = "Component", y = "Variance explained (%)",
                  title = "Scree plot",
                  subtitle = "Bars: per-PC  |  Line: cumulative") +
    theme_libs()
}

#' Plot PLS-DA Results
#'
#' @param plsda A `libs_plsda` object from [ls_plsda()].
#' @param type Character. `"scores"` (LV1 vs LV2), `"confusion"` (heatmap),
#'   or `"vip"` (VIP scores spectrum). Default `"scores"`.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' \donttest{
#' if (requireNamespace("pls", quietly = TRUE)) {
#'   plsda <- ls_plsda(ds, "tissue", n_components = 3, validation = "none")
#'   ls_plot_plsda(plsda, type = "scores")
#' }
#' }
#' @export
ls_plot_plsda <- function(plsda, type = "scores") {
  if (!inherits(plsda, "libs_plsda")) {
    cli::cli_abort("{.arg plsda} must be a {.cls libs_plsda} object.")
  }
  type <- match.arg(type, c("scores", "confusion", "vip"))

  if (type == "scores") {
    scores <- plsda$model$scores[, 1:min(2, ncol(plsda$model$scores)),
                                 drop = FALSE]
    df <- data.frame(
      LV1 = scores[, 1],
      LV2 = if (ncol(scores) >= 2) scores[, 2] else rep(0, nrow(scores)),
      class = plsda$predictions$observed
    )
    return(ggplot2::ggplot(df, ggplot2::aes(x = .data$LV1, y = .data$LV2,
                                            colour = .data$class)) +
             ggplot2::geom_point(size = 2.5, alpha = 0.8) +
             ggplot2::stat_ellipse(level = 0.95, alpha = 0.5) +
             ggplot2::scale_colour_brewer(palette = "Set1") +
             ggplot2::labs(x = "LV1", y = "LV2",
                           title = sprintf("PLS-DA Scores (accuracy %.1f%%)",
                                           plsda$accuracy * 100)) +
             theme_libs())
  }
  if (type == "confusion") {
    cm <- plsda$confusion_matrix
    df <- as.data.frame(cm)
    names(df) <- c("observed", "predicted", "count")
    return(ggplot2::ggplot(df, ggplot2::aes(x = .data$predicted,
                                            y = .data$observed,
                                            fill = .data$count)) +
             ggplot2::geom_tile() +
             ggplot2::geom_text(ggplot2::aes(label = .data$count),
                                colour = "white") +
             ggplot2::scale_fill_viridis_c(option = "viridis") +
             ggplot2::labs(x = "Predicted", y = "Observed",
                           title = "Confusion matrix") +
             theme_libs())
  }
  # VIP
  vip <- .vip_scores(plsda$model, plsda$n_components)
  wl <- plsda$model$scale
  df <- data.frame(wavelength_idx = seq_along(vip), vip = vip)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$wavelength_idx,
                                   y = .data$vip)) +
    ggplot2::geom_line(colour = "#2E86AB") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "red") +
    ggplot2::labs(x = "Wavelength channel", y = "VIP score",
                  title = "Variable Importance in Projection") +
    theme_libs()
}

# VIP scores for PLS models
#' @keywords internal
#' @noRd
.vip_scores <- function(object, ncomp) {
  if (!requireNamespace("pls", quietly = TRUE)) {
    return(rep(NA_real_, length(object$coefficients)))
  }
  W <- object$loading.weights
  Q <- object$Yloadings
  TT <- object$scores
  p <- nrow(W)
  h <- ncomp
  q <- ncol(Q)
  # SS per component: (Q[,,h]^2) * (TT[,h]'*TT[,h])
  SS <- rep(0, h)
  for (k in seq_len(h)) {
    SS[k] <- sum(Q[, k]^2) * sum(TT[, k]^2)
  }
  vip <- numeric(p)
  for (i in seq_len(p)) {
    s <- 0
    for (k in seq_len(h)) {
      w_ik <- W[i, k]
      denom <- sum(W[, k]^2)
      if (denom == 0) next
      s <- s + SS[k] * (w_ik^2 / denom)
    }
    vip[i] <- sqrt(p * s / sum(SS))
  }
  vip
}

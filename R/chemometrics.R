# Chemometrics: PCA, PLS-DA, clustering, classifiers
# -----------------------------------------------------------------------------

#' Principal Component Analysis of LIBS Spectra
#'
#' Performs PCA on the intensity matrix of a [ls_dataset()].
#'
#' @param dataset A [ls_dataset()] object.
#' @param n_components Integer. Number of PCs to retain. Default 10.
#' @param scale Logical. Scale variables to unit variance. Default `TRUE`.
#' @param center Logical. Center variables. Default `TRUE`.
#'
#' @return A list with S3 class `libs_pca` containing:
#'   `scores` (matrix: spectra x PCs), `loadings` (matrix: channels x PCs),
#'   `sdev`, `variance_explained` (fraction), `cumulative_variance`,
#'   `wavelength`, `sample_info`, `pca_obj`.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' pca <- ls_pca(ds, n_components = 5)
#' pca$variance_explained
#'
#' @export
ls_pca <- function(dataset, n_components = 10, scale = TRUE, center = TRUE) {
  .validate_dataset(dataset)
  X <- dataset$intensity_matrix
  # Remove zero-variance columns if scale=TRUE to avoid errors
  if (isTRUE(scale)) {
    col_sds <- apply(X, 2, stats::sd)
    nonzero <- col_sds > 0
    if (!all(nonzero)) {
      X <- X[, nonzero, drop = FALSE]
      wl_used <- dataset$wavelength[nonzero]
    } else {
      wl_used <- dataset$wavelength
    }
  } else {
    wl_used <- dataset$wavelength
  }
  pca <- stats::prcomp(X, center = center, scale. = scale)
  n_components <- min(n_components, ncol(pca$x))
  var_expl <- pca$sdev^2 / sum(pca$sdev^2)

  # Pad loadings back to full wavelength axis
  full_loadings <- matrix(0,
                          nrow = dataset$n_channels,
                          ncol = n_components)
  if (length(wl_used) != dataset$n_channels) {
    idx_full <- match(wl_used, dataset$wavelength)
    full_loadings[idx_full, ] <- pca$rotation[, seq_len(n_components), drop = FALSE]
  } else {
    full_loadings <- pca$rotation[, seq_len(n_components), drop = FALSE]
  }

  structure(
    list(
      scores = pca$x[, seq_len(n_components), drop = FALSE],
      loadings = full_loadings,
      sdev = pca$sdev[seq_len(n_components)],
      variance_explained = var_expl[seq_len(n_components)],
      cumulative_variance = cumsum(var_expl)[seq_len(n_components)],
      wavelength = dataset$wavelength,
      sample_info = dataset$sample_info,
      n_components = n_components,
      pca_obj = pca
    ),
    class = "libs_pca"
  )
}

#' @export
print.libs_pca <- function(x, ...) {
  cli::cli_inform(c(
    "{.cls libs_pca}",
    "*" = "Components: {x$n_components}",
    "*" = "Spectra: {nrow(x$scores)}",
    "*" = "Variance (first 3): {paste(round(x$variance_explained[seq_len(min(3, x$n_components))] * 100, 1), collapse = ', ')}%",
    "*" = "Cumulative (first 3): {round(x$cumulative_variance[min(3, x$n_components)] * 100, 1)}%"
  ))
  invisible(x)
}

#' PLS-DA Classification
#'
#' Partial Least Squares Discriminant Analysis for class prediction.
#'
#' @param dataset A [ls_dataset()] object.
#' @param grouping Character. Column name in `sample_info` with class labels.
#' @param n_components Integer. Number of PLS latent variables. Default 5.
#' @param validation Character. `"CV"` (10-fold), `"LOO"`, or `"none"`.
#'   Default `"CV"`.
#'
#' @return An S3 object of class `libs_plsda` with elements `model`,
#'   `predictions` (tibble with predicted/observed class), `confusion_matrix`,
#'   `accuracy`, `class_labels`, `n_components`, `sample_info`.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' \donttest{
#' if (requireNamespace("pls", quietly = TRUE)) {
#'   plsda <- ls_plsda(ds, "tissue", n_components = 3, validation = "none")
#'   plsda$accuracy
#' }
#' }
#' @export
ls_plsda <- function(dataset, grouping, n_components = 5, validation = "CV") {
  .validate_dataset(dataset)
  .require_pkg("pls", "PLS-DA requires the {.pkg pls} package.")
  validation <- match.arg(validation, c("CV", "LOO", "none"))

  if (!grouping %in% names(dataset$sample_info)) {
    cli::cli_abort("{.arg grouping} column {.val {grouping}} not found in sample_info.")
  }
  labels <- as.factor(dataset$sample_info[[grouping]])
  class_labels <- levels(labels)
  k <- length(class_labels)
  if (k < 2) {
    cli::cli_abort("Need at least 2 classes; got {k}.")
  }

  Y <- stats::model.matrix(~ labels - 1)
  colnames(Y) <- class_labels
  X <- dataset$intensity_matrix
  n_components <- min(n_components, nrow(X) - 1, ncol(X))

  val_arg <- switch(validation, CV = "CV", LOO = "LOO", none = "none")

  df <- data.frame(Y = I(Y))
  df$X <- X
  fit <- pls::plsr(Y ~ X, ncomp = n_components, data = df,
                   validation = val_arg,
                   segments = if (validation == "CV") min(10, nrow(X)) else NULL)
  pred_mat <- stats::predict(fit, ncomp = n_components)
  pred_mat <- pred_mat[, , 1]
  pred_class <- class_labels[apply(pred_mat, 1, which.max)]
  obs_class <- as.character(labels)
  cm <- table(observed = obs_class, predicted = pred_class)
  acc <- sum(obs_class == pred_class) / length(obs_class)

  structure(
    list(
      model = fit,
      predictions = tibble::tibble(
        sample_id = dataset$sample_info$sample_id,
        observed = obs_class,
        predicted = pred_class
      ),
      confusion_matrix = cm,
      accuracy = acc,
      class_labels = class_labels,
      n_components = n_components,
      sample_info = dataset$sample_info
    ),
    class = "libs_plsda"
  )
}

#' @export
print.libs_plsda <- function(x, ...) {
  cli::cli_inform(c(
    "{.cls libs_plsda}",
    "*" = "Classes ({length(x$class_labels)}): {paste(x$class_labels, collapse = ', ')}",
    "*" = "Components: {x$n_components}",
    "*" = "Accuracy: {round(x$accuracy * 100, 1)}%"
  ))
  invisible(x)
}

#' Cluster LIBS Spectra
#'
#' Unsupervised clustering via k-means, hierarchical, or DBSCAN (DBSCAN
#' requires the `dbscan` package in Suggests).
#'
#' @param dataset A [ls_dataset()] object.
#' @param method Character. `"kmeans"`, `"hclust"`, or `"dbscan"`.
#'   Default `"kmeans"`.
#' @param k Integer. Number of clusters for kmeans/hclust. Default 3.
#' @param ... Additional arguments passed to the underlying clustering
#'   function.
#'
#' @return An S3 object of class `libs_clusters` with elements `cluster`,
#'   `method`, `k`, `silhouette` (mean silhouette score when computable),
#'   `sample_info`.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' cl <- ls_cluster(ds, method = "kmeans", k = 5)
#' cl$cluster
#' @export
ls_cluster <- function(dataset, method = "kmeans", k = 3, ...) {
  .validate_dataset(dataset)
  method <- match.arg(method, c("kmeans", "hclust", "dbscan"))
  X <- dataset$intensity_matrix
  # Standardize
  X_scaled <- scale(X)
  X_scaled[is.na(X_scaled)] <- 0

  if (method == "kmeans") {
    km <- stats::kmeans(X_scaled, centers = k, nstart = 10, ...)
    clusters <- km$cluster
    centers <- km$centers
  } else if (method == "hclust") {
    d <- stats::dist(X_scaled)
    hc <- stats::hclust(d, method = "ward.D2")
    clusters <- stats::cutree(hc, k = k)
    centers <- NULL
  } else {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      cli::cli_abort("DBSCAN requires the {.pkg dbscan} package.")
    }
    db <- dbscan::dbscan(X_scaled, ...)
    clusters <- db$cluster
    k <- length(unique(clusters[clusters > 0]))
    centers <- NULL
  }

  sil <- NA_real_
  if (k >= 2 && length(unique(clusters)) >= 2) {
    sil <- tryCatch({
      d <- stats::dist(X_scaled)
      .silhouette_score(d, clusters)
    }, error = function(e) NA_real_)
  }

  structure(
    list(
      cluster = clusters,
      method = method,
      k = k,
      centers = centers,
      silhouette = sil,
      sample_info = dataset$sample_info
    ),
    class = "libs_clusters"
  )
}

#' @export
print.libs_clusters <- function(x, ...) {
  cli::cli_inform(c(
    "{.cls libs_clusters}",
    "*" = "Method: {.val {x$method}}",
    "*" = "k: {x$k}",
    "*" = "Silhouette: {if (is.na(x$silhouette)) 'NA' else round(x$silhouette, 3)}"
  ))
  invisible(x)
}

#' Train a Classifier
#'
#' Trains SVM or Random Forest on spectral intensity features.
#'
#' @param dataset A [ls_dataset()] with labelled samples.
#' @param grouping Character. Column name in `sample_info` with class labels.
#' @param method Character. `"svm"` or `"rf"`. Default `"svm"`.
#' @param ... Additional arguments passed to the underlying trainer.
#'
#' @return An S3 object of class `libs_classifier`.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' \donttest{
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   clf <- ls_train_classifier(ds, "tissue", method = "svm")
#'   clf$accuracy
#' }
#' }
#' @export
ls_train_classifier <- function(dataset, grouping, method = "svm", ...) {
  .validate_dataset(dataset)
  method <- match.arg(method, c("svm", "rf"))
  if (!grouping %in% names(dataset$sample_info)) {
    cli::cli_abort("{.arg grouping} column not found in sample_info.")
  }
  labels <- as.factor(dataset$sample_info[[grouping]])
  X <- dataset$intensity_matrix

  if (method == "svm") {
    .require_pkg("e1071", "SVM training requires {.pkg e1071}.")
    df <- data.frame(X)
    df$.class <- labels
    fit <- e1071::svm(.class ~ ., data = df, probability = TRUE, ...)
    pred <- stats::predict(fit, df)
    acc <- mean(pred == labels)
  } else {
    .require_pkg("ranger", "Random Forest requires {.pkg ranger}.")
    df <- data.frame(X)
    df$.class <- labels
    fit <- ranger::ranger(.class ~ ., data = df,
                          classification = TRUE,
                          probability = TRUE, ...)
    pred_probs <- fit$predictions
    pred <- apply(pred_probs, 1, function(r) levels(labels)[which.max(r)])
    acc <- mean(pred == labels)
  }

  cm <- table(observed = labels, predicted = pred)

  structure(
    list(
      model = fit,
      method = method,
      accuracy = acc,
      confusion_matrix = cm,
      class_labels = levels(labels),
      sample_info = dataset$sample_info
    ),
    class = "libs_classifier"
  )
}

#' @export
print.libs_classifier <- function(x, ...) {
  cli::cli_inform(c(
    "{.cls libs_classifier}",
    "*" = "Method: {.val {x$method}}",
    "*" = "Classes ({length(x$class_labels)}): {paste(x$class_labels, collapse = ', ')}",
    "*" = "Training accuracy: {round(x$accuracy * 100, 1)}%"
  ))
  invisible(x)
}

#' Classify Unknown Spectra
#'
#' Applies a trained classifier (from [ls_plsda()] or [ls_train_classifier()])
#' to new spectra.
#'
#' @param model A `libs_plsda` or `libs_classifier` object.
#' @param new_data A [ls_dataset()] of unknown spectra.
#'
#' @return A [tibble::tibble()] with columns `sample_id`, `predicted_class`,
#'   and optionally `probability`.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' \donttest{
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   clf <- ls_train_classifier(ds[1:30], "tissue", method = "svm")
#'   ls_classify(clf, ds[31:50])
#' }
#' }
#' @export
ls_classify <- function(model, new_data) {
  .validate_dataset(new_data)
  X <- new_data$intensity_matrix
  sids <- new_data$sample_info$sample_id

  if (inherits(model, "libs_plsda")) {
    pred_mat <- stats::predict(model$model, newdata = X, ncomp = model$n_components)
    pred_mat <- pred_mat[, , 1]
    pred <- model$class_labels[apply(pred_mat, 1, which.max)]
    probs <- apply(pred_mat, 1, max)
    return(tibble::tibble(
      sample_id = sids,
      predicted_class = pred,
      probability = probs
    ))
  }
  if (inherits(model, "libs_classifier")) {
    if (model$method == "svm") {
      df <- data.frame(X)
      pred <- stats::predict(model$model, df, probability = TRUE)
      probs <- attr(pred, "probabilities")
      prob_max <- if (is.null(probs)) rep(NA_real_, length(pred)) else apply(probs, 1, max)
      return(tibble::tibble(
        sample_id = sids,
        predicted_class = as.character(pred),
        probability = prob_max
      ))
    }
    if (model$method == "rf") {
      df <- data.frame(X)
      prd <- stats::predict(model$model, data = df)
      probs <- prd$predictions
      pred <- model$class_labels[apply(probs, 1, which.max)]
      prob_max <- apply(probs, 1, max)
      return(tibble::tibble(
        sample_id = sids,
        predicted_class = pred,
        probability = prob_max
      ))
    }
  }
  cli::cli_abort("Unsupported model class: {.val {class(model)[1]}}")
}

# Silhouette score (mean across samples)
#' @keywords internal
#' @noRd
.silhouette_score <- function(d, clusters) {
  d <- as.matrix(d)
  n <- length(clusters)
  if (n < 3) return(NA_real_)
  scores <- numeric(n)
  for (i in seq_len(n)) {
    ci <- clusters[i]
    same <- which(clusters == ci & seq_len(n) != i)
    a_i <- if (length(same) > 0) mean(d[i, same]) else 0
    other_clusters <- setdiff(unique(clusters), ci)
    if (length(other_clusters) == 0) {
      scores[i] <- 0
      next
    }
    b_vec <- vapply(other_clusters, function(oc) {
      mem <- which(clusters == oc)
      if (length(mem) == 0) Inf else mean(d[i, mem])
    }, numeric(1))
    b_i <- min(b_vec)
    scores[i] <- (b_i - a_i) / max(a_i, b_i, 1e-12)
  }
  mean(scores, na.rm = TRUE)
}

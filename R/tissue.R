# Tissue-specific classification and discrimination
# -----------------------------------------------------------------------------

# Canonical lines for tissue analysis
.tissue_lines <- c(
  Ca = 393.37, Na = 589.00, K = 766.49, Mg = 279.55,
  Fe = 371.99, Zn = 213.86, Cu = 324.75, P = 213.62,
  C = 247.86, H = 656.28
)

# Rule-based tissue discrimination using elemental ratios.
# Returns best-matching tissue label and a confidence score.
#' @keywords internal
#' @noRd
.classify_by_ratio <- function(s) {
  areas <- vapply(names(.tissue_lines), function(el) {
    wn <- .tissue_lines[[el]]
    if (wn < min(s$wavelength) || wn > max(s$wavelength)) return(0)
    ls_peak_area(s, wn, 1)
  }, numeric(1))
  names(areas) <- names(.tissue_lines)
  # Guard against zeros
  safe <- pmax(areas, 1e-6)

  # Rules (relative)
  ca_score <- safe["Ca"] / (safe["Na"] + safe["K"] + safe["Fe"])
  fe_score <- safe["Fe"] / (safe["Na"] + safe["K"] + safe["Ca"])
  k_score <- safe["K"] / (safe["Na"] + safe["Ca"])
  zn_score <- safe["Zn"] / (safe["Na"] + safe["K"])
  c_score <- safe["C"] / (safe["Ca"] + safe["Fe"] + safe["K"])

  scores <- c(
    bone = as.numeric(ca_score),
    liver = as.numeric(fe_score),
    kidney = as.numeric(k_score * (safe["Na"] / safe["Ca"])),
    muscle = as.numeric(k_score),
    fat = as.numeric(c_score)
  )
  scores <- scores / sum(abs(scores) + 1e-12)
  best <- names(scores)[which.max(scores)]
  list(
    tissue = best,
    confidence = max(scores),
    scores = scores,
    areas = areas
  )
}

#' Tissue Classification via LIBS
#'
#' Classifies tissue types based on LIBS spectral signatures. Implements
#' three methods:
#' * `"ratio"` — rule-based using canonical elemental ratios (Ca/Na, Fe/Na,
#'   Zn/Na, K/Na, C/Ca) derived from tissue biochemistry. No training data
#'   required.
#' * `"plsda"` — fits a PLS-DA model on the supplied labelled `reference`
#'   dataset and predicts class for `x`.
#' * `"svm"` — fits an SVM on `reference` and predicts class.
#'
#' @param x A [ls_spectrum()] or [ls_dataset()] object (unknowns).
#' @param method Character. `"ratio"`, `"plsda"`, or `"svm"`. Default `"ratio"`.
#' @param reference A [ls_dataset()] with labelled tissue reference spectra
#'   (needed for `"plsda"` and `"svm"`).
#' @param group_col Character. Column in reference `sample_info` with tissue
#'   labels. Default `"tissue"`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A [tibble::tibble()] with columns `sample_id`, `predicted_tissue`,
#'   `confidence`.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' ls_tissue_classify(ds[1:3])
#' @export
ls_tissue_classify <- function(x, method = "ratio", reference = NULL,
                               group_col = "tissue", verbose = TRUE) {
  method <- match.arg(method, c("ratio", "plsda", "svm"))

  if (inherits(x, "libs_spectrum")) {
    if (method == "ratio") {
      r <- .classify_by_ratio(x)
      sid <- as.character(.meta_get(x$metadata, "sample_id", "unknown"))
      return(tibble::tibble(
        sample_id = sid,
        predicted_tissue = r$tissue,
        confidence = r$confidence
      ))
    }
    # Wrap single spectrum as dataset
    x_ds <- ls_dataset(list(x))
  } else {
    .validate_dataset(x)
    x_ds <- x
  }

  if (method == "ratio") {
    rows <- lapply(x_ds$spectra, function(s) {
      r <- .classify_by_ratio(s)
      sid <- as.character(.meta_get(s$metadata, "sample_id", NA))
      data.frame(sample_id = sid, predicted_tissue = r$tissue,
                 confidence = r$confidence, stringsAsFactors = FALSE)
    })
    return(tibble::as_tibble(do.call(rbind, rows)))
  }

  if (is.null(reference)) {
    cli::cli_abort("{.arg reference} is required for method {.val {method}}.")
  }
  .validate_dataset(reference)
  if (!group_col %in% names(reference$sample_info)) {
    cli::cli_abort("Column {.val {group_col}} missing from {.arg reference$sample_info}.")
  }

  if (method == "plsda") {
    plsda <- ls_plsda(reference, grouping = group_col,
                      n_components = 5, validation = "none")
    return(ls_classify(plsda, x_ds) |>
             dplyr::rename(predicted_tissue = "predicted_class",
                           confidence = "probability"))
  }
  if (method == "svm") {
    clf <- ls_train_classifier(reference, grouping = group_col, method = "svm")
    return(ls_classify(clf, x_ds) |>
             dplyr::rename(predicted_tissue = "predicted_class",
                           confidence = "probability"))
  }
}

#' Tissue Discrimination Analysis
#'
#' Identifies emission channels that best discriminate two tissue types via
#' per-wavelength t-tests with Benjamini-Hochberg FDR correction, plus log2
#' fold-change. Optionally matches significant channels to NIST elements.
#'
#' @param dataset A [ls_dataset()] with tissue labels.
#' @param group_col Character. Label column in `sample_info`.
#' @param group_a Character. First group label.
#' @param group_b Character. Second group label.
#' @param method Character. `"t_test"` (default) or `"fold_change"`.
#' @param tolerance_nm Numeric. Matching tolerance when annotating elements.
#'   Default 0.2.
#'
#' @return A [tibble::tibble()] with columns `wavelength_nm`, `mean_a`,
#'   `mean_b`, `p_value`, `fold_change`, `fdr`, `significant`, `element`.
#'
#' @examples
#' ds <- ls_example_data("tissue")
#' head(ls_tissue_discriminate(ds, "tissue", "bone", "muscle"))
#' @export
ls_tissue_discriminate <- function(dataset, group_col, group_a, group_b,
                                   method = "t_test", tolerance_nm = 0.2) {
  .validate_dataset(dataset)
  method <- match.arg(method, c("t_test", "fold_change"))
  if (!group_col %in% names(dataset$sample_info)) {
    cli::cli_abort("{.arg group_col} not found in sample_info.")
  }
  labs <- as.character(dataset$sample_info[[group_col]])
  ia <- which(labs == group_a)
  ib <- which(labs == group_b)
  if (length(ia) < 2 || length(ib) < 2) {
    cli::cli_abort("Each group must have at least 2 samples.")
  }
  Xa <- dataset$intensity_matrix[ia, , drop = FALSE]
  Xb <- dataset$intensity_matrix[ib, , drop = FALSE]
  mean_a <- colMeans(Xa)
  mean_b <- colMeans(Xb)
  fc <- log2((mean_a + 1) / (mean_b + 1))

  if (method == "t_test") {
    p_vals <- vapply(seq_len(ncol(Xa)), function(j) {
      tryCatch(stats::t.test(Xa[, j], Xb[, j])$p.value,
               error = function(e) NA_real_)
    }, numeric(1))
  } else {
    p_vals <- rep(NA_real_, ncol(Xa))
  }
  fdr <- if (all(is.na(p_vals))) p_vals else stats::p.adjust(p_vals, method = "BH")
  sig <- ifelse(is.na(fdr), FALSE, fdr < 0.05 & abs(fc) > 1)

  # Annotate elements
  el_db <- ls_element_db()
  elements <- rep(NA_character_, length(dataset$wavelength))
  for (i in seq_along(dataset$wavelength)) {
    dev <- dataset$wavelength[i] - el_db$wavelength_nm
    within <- abs(dev) <= tolerance_nm
    if (!any(within)) next
    cand <- el_db[within, ]
    score <- (1 - abs(dev[within]) / tolerance_nm) * cand$aki
    elements[i] <- cand$element[which.max(score)]
  }

  out <- tibble::tibble(
    wavelength_nm = dataset$wavelength,
    mean_a = mean_a,
    mean_b = mean_b,
    p_value = p_vals,
    fold_change = fc,
    fdr = fdr,
    significant = sig,
    element = elements
  )
  out[order(out$fdr, -abs(out$fold_change)), ]
}

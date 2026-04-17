# Declare variables used in NSE (tidy eval) to silence R CMD check
utils::globalVariables(c(
  "wavelength_nm", "intensity", "element", "sample_id",
  "concentration", "group", "x_pos", "y_pos", "value",
  "ionization", "aki", "ei_ev", "ek_ev", "persistent",
  "PC1", "PC2", "LV1", "LV2", "predicted", "observed",
  "fitted_value", "residual", "snr", "prominence",
  "fwhm_nm", "area", "wavelength", "deviation_nm",
  "confidence", "nist_wavelength_nm", "p_value",
  "fold_change", "fdr", "significant", "cluster",
  "loading", "score", "component", "variance_explained",
  "term", "configuration", "class_predicted", "class_true",
  "count", "prob", "tissue", "predicted_tissue",
  "shot", "rep", "spectrum_id", "idx"
))

test_that("ls_quantify returns tibble with expected columns", {
  ds <- ls_example_data("calibration", n_channels = 256)
  conc <- ds$sample_info$concentration
  cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
  out <- ls_quantify(cal, ds)
  expect_s3_class(out, "tbl_df")
  expect_true(all(c("sample_id", "element", "concentration",
                    "unit", "below_lod", "below_loq") %in% names(out)))
  expect_equal(nrow(out), ds$n_spectra)
})

test_that("ls_quantify on single spectrum works", {
  ds <- ls_example_data("calibration", n_channels = 256)
  conc <- ds$sample_info$concentration
  cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
  out <- ls_quantify(cal, ds$spectra[[1]])
  expect_equal(nrow(out), 1)
})

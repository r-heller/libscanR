test_that("ls_calibrate produces linear fit with high R^2", {
  ds <- ls_example_data("calibration", n_channels = 512)
  conc <- ds$sample_info$concentration
  cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
  expect_s3_class(cal, "libs_calibration")
  expect_gt(cal$r_squared, 0.8)
  expect_true(!is.null(cal$lod))
})

test_that("predict.libs_calibration returns numeric predictions", {
  ds <- ls_example_data("calibration", n_channels = 256)
  conc <- ds$sample_info$concentration
  cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
  pred <- predict(cal, newdata = cal$intensities)
  expect_type(pred, "double")
  expect_length(pred, length(cal$intensities))
})

test_that("ls_lod and ls_loq return positive values", {
  ds <- ls_example_data("calibration", n_channels = 256)
  conc <- ds$sample_info$concentration
  cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
  lod <- ls_lod(cal)
  loq <- ls_loq(cal)
  expect_gt(lod, 0)
  expect_gt(loq, lod)
})

test_that("ls_saha_boltzmann returns a tibble with temperature", {
  skip_on_cran()
  s <- ls_simulate_spectrum(
    elements = c(Ca = 10000, Fe = 500, Na = 2000),
    seed = 3, n_channels = 1024
  )
  s <- ls_baseline(s, method = "snip", iterations = 30)
  res <- ls_saha_boltzmann(
    s, elements = c("Ca", "Fe"),
    lines_nm = list(Ca = c(422.673, 445.478, 487.813),
                    Fe = c(371.994, 404.581, 438.354)),
    verbose = FALSE
  )
  expect_s3_class(res, "tbl_df")
  expect_true("temperature_k" %in% names(res))
})

test_that("print.libs_calibration does not error", {
  ds <- ls_example_data("calibration", n_channels = 256)
  conc <- ds$sample_info$concentration
  cal <- ls_calibrate(ds, "Ca", 393.37, conc, verbose = FALSE)
  expect_message(print(cal), "libs_calibration")
})

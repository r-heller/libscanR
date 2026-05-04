test_that("ls_build_map constructs libs_map", {
  ds <- ls_example_data("spatial", n_channels = 128)
  m <- ls_build_map(ds, "Ca", 393.37)
  expect_s3_class(m, "libs_map")
  expect_equal(length(m$values), ds$n_spectra)
  expect_equal(m$element, "Ca")
})

test_that("ls_build_map builds a grid for regular coordinates", {
  ds <- ls_example_data("spatial", n_channels = 128)
  m <- ls_build_map(ds, "Ca", 393.37)
  expect_false(is.null(m$grid))
  expect_equal(dim(m$grid), c(20, 20))
})

test_that("ls_map_elements returns list of maps", {
  ds <- ls_example_data("spatial", n_channels = 128)
  ms <- ls_map_elements(ds, c("Ca", "Fe"),
                        c(Ca = 393.37, Fe = 371.99))
  expect_type(ms, "list")
  expect_length(ms, 2)
  expect_s3_class(ms[["Ca"]], "libs_map")
})

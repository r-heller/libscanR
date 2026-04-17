test_that("ls_run_app function exists and app files are installed", {
  expect_true(is.function(ls_run_app))
  # When run via devtools::load_all, system.file returns "" — skip
  # the path check in that case.
  path <- system.file("shiny", "libscanR", package = "libscanR")
  if (nzchar(path)) {
    expect_true(file.exists(file.path(path, "app.R")))
  }
})

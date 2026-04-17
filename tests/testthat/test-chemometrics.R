test_that("ls_pca returns libs_pca object", {
  ds <- ls_example_data("tissue", n_channels = 128)
  pca <- ls_pca(ds, n_components = 4)
  expect_s3_class(pca, "libs_pca")
  expect_equal(pca$n_components, 4)
  expect_length(pca$variance_explained, 4)
  expect_true(all(pca$cumulative_variance >= pca$variance_explained))
})

test_that("ls_plsda runs with pls package", {
  skip_if_not_installed("pls")
  ds <- ls_example_data("tissue", n_channels = 128)
  plsda <- ls_plsda(ds, grouping = "tissue",
                    n_components = 3, validation = "none")
  expect_s3_class(plsda, "libs_plsda")
  expect_gte(plsda$accuracy, 0)
  expect_lte(plsda$accuracy, 1)
})

test_that("ls_cluster (kmeans) assigns cluster labels", {
  ds <- ls_example_data("tissue", n_channels = 128)
  cl <- ls_cluster(ds, method = "kmeans", k = 5)
  expect_s3_class(cl, "libs_clusters")
  expect_equal(length(cl$cluster), ds$n_spectra)
  expect_lte(length(unique(cl$cluster)), 5)
})

test_that("ls_cluster (hclust) works", {
  ds <- ls_example_data("tissue", n_channels = 128)[1:20]
  cl <- ls_cluster(ds, method = "hclust", k = 3)
  expect_s3_class(cl, "libs_clusters")
})

test_that("ls_train_classifier (svm) works", {
  skip_if_not_installed("e1071")
  ds <- ls_example_data("tissue", n_channels = 64)
  clf <- ls_train_classifier(ds, "tissue", method = "svm")
  expect_s3_class(clf, "libs_classifier")
  expect_gte(clf$accuracy, 0)
})

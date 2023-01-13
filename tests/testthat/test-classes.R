test_that("All modules are the correct S4 classes", {
  m <- create_model_Scaled(10, 3, 1:10, create_Lmatrix(1:10, 3, 10))
  expect_s4_class(m, "Rcpp_Scaled")
  m <- create_model_Unscaled(10, 3, 1:10, create_Lmatrix(1:10, 3, 10))
  expect_s4_class(m, "Rcpp_Unscaled")
  m <- create_model_Unscaled_nuts(10, 3, 2, 1:10, create_Lmatrix(1:10, 3, 10))
  expect_s4_class(m, "Rcpp_Unscaled_nuts")
})

test_that("All classes input is correct", {
  set.seed(1234)
  nb_s <- sample(seq(2, 1e1, by = 2), 1)
  fw <- create_niche_model(nb_s, .3)
  nb_n <- sample(seq(2, 1e1), 1)
  # scaled
  m <- create_model_Scaled(nb_s, nb_s / 2, 1:nb_s, fw)
  expect_equal(m$nb_s, nb_s)
  expect_equal(m$nb_b, nb_s / 2)
  # unscaled
  m <- create_model_Unscaled(nb_s, nb_s / 2, 1:nb_s, fw)
  expect_equal(m$nb_s, nb_s)
  expect_equal(m$nb_b, nb_s / 2)
  # nutrients
  m <- create_model_Unscaled_nuts(nb_s, nb_s / 2, nb_n, 1:nb_s, fw)
  expect_equal(m$nb_s, nb_s)
  expect_equal(m$nb_b, nb_s / 2)
  expect_equal(m$nb_n, nb_n)
})

test_that("Initializations are class specific", {
  set.seed(1234)
  nb_s <- sample(seq(2, 1e1, by = 2), 1)
  fw <- create_niche_model(nb_s, .3)
  nb_n <- sample(seq(2, 1e1), 1)
  # scaled
  m <- create_model_Scaled(nb_s, nb_s / 2, 1:nb_s, fw)
  expect_error(initialise_default_Unscaled(m))
  expect_error(initialise_default_Unscaled_nuts(m))
  # unscaled
  m <- create_model_Unscaled(nb_s, nb_s / 2, 1:nb_s, fw)
  expect_error(initialise_default_Scaled(m))
  expect_error(initialise_default_Unscaled_nuts(m))
  # nutrients
  m <- create_model_Unscaled_nuts(nb_s, nb_s / 2, nb_n, 1:nb_s, fw)
  expect_error(initialise_default_Scaled(m))
  expect_error(initialise_default_Unscaled(m))
})

test_that("Initializations work properly", {
  set.seed(1234)
  nb_s <- sample(seq(2, 1e1, by = 2), 1)
  fw <- create_niche_model(nb_s, .3)
  nb_n <- sample(seq(2, 1e1), 1)
  # scaled
  m <- create_model_Scaled(nb_s, nb_s / 2, 1:nb_s, fw)
  expect_equal(sum(m$alpha), 0)
  expect_equal(sum(m$X), 0)
  m <- initialise_default_Scaled(m)
  expect_gt(sum(m$alpha), 0)
  expect_gt(sum(m$X), 0)
  # unscaled
  m <- create_model_Unscaled(nb_s, nb_s / 2, 1:nb_s, fw)
  expect_equal(sum(m$alpha), 0)
  expect_equal(sum(m$X), 0)
  m <- initialise_default_Unscaled(m)
  expect_gt(sum(m$alpha), 0)
  expect_gt(sum(m$X), 0)
  # nutrients
  m <- create_model_Unscaled_nuts(nb_s, nb_s / 2, nb_n, 1:nb_s, fw)
  expect_equal(sum(m$D), 0)
  expect_equal(sum(m$X), 0)
  expect_equal(sum(m$V), 0)
  expect_equal(sum(m$S), 0)
  m <- initialise_default_Unscaled_nuts(m, L.mat = fw)
  expect_gt(sum(m$D), 0)
  expect_gt(sum(m$X), 0)
  expect_gt(sum(m$V), 0)
  expect_gt(sum(m$S), 0)
})

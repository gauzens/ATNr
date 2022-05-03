test_that("sorting matrix works properly", {
  set.seed(123)
  fw <- matrix(as.numeric(runif(100) > .9), 10, 10)
  bm <- runif(10)
  sp <- letters[1:10]
  names(bm) <- sp
  rownames(fw) <- sp
  colnames(fw) <- sp
  sorted <- sort_input(bm, fw)
  # check in degrees
  expect_identical(colSums(sorted$food.web)[order(colnames(sorted$food.web))],
                   colSums(fw))
  # check body masses
  expect_identical(sorted$body.mass[order(names(sorted$body.mass))],
                   bm)
  # check total number of links
  expect_equal(sum(fw), sum(sorted$food.web))
})


test_that("covarianceSOV", {
  expect_error(covarianceSOV(diag(1,2), c(0,0), c(1, -1)));
  expect_error(covarianceSOV(diag(1,2), c(0,0,0), c(1,1)));
  expect_lte(covarianceSOV(diag(1,2), c(0,0), c(1,1)), 1);
  expect_gte(covarianceSOV(diag(1,2), c(0,0), c(1,1)), 0);
  expect_equal(covarianceSOV(diag(1,2), c(0,0), c(1,1)), precisionSOV(diag(1,2), c(0,0), c(1,1)))
  expect_equal(covarianceSOV(diag(7,100), rep(0,100), rep(1,100)), precisionSOV(diag(1/7,100), rep(0,100), rep(1,100)))
})

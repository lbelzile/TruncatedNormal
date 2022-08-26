test_that("precisionSOV", {
  expect_error(precisionSOV(diag(1,2), c(0,0), c(1, -1)));
  expect_error(precisionSOV(diag(1,2), c(0,0,0), c(1,1)));
  expect_lte(precisionSOV(diag(1,2), c(0,0), c(1,1)), 1);
  expect_gte(precisionSOV(diag(1,2), c(0,0), c(1,1)), 0);
  expect_equal(precisionSOV(diag(1,2), c(0,0), c(1,1)), covarianceSOV(diag(1,2), c(0,0), c(1,1)))
  expect_equal(precisionSOV(diag(7,100), rep(0,100), rep(1,100)), covarianceSOV(diag(1/7,100), rep(0,100), rep(1,100)))
  
})

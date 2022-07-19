test_that("trnd2", {
  expect_error(trnd2(c(-1,-1), c(0,0,0)));
  expect_error(trnd2(c(-1,-1), c(0, -2)));
  expect_vector(trnd2(c(-1,0), c(0,1)));
  expect_lte(trnd2(-1,0), 0);
  expect_gte(trnd2(-1,0), -1);
})

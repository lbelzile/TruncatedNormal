test_that("tn", {
  expect_vector(tn(c(0,1,2,3), c(1,2,3,4)));
  expect_length(tn(c(0,1), c(1,2)), 2);
  expect_lte(tn(0,1), 1);
  expect_gte(tn(0,1), 0);
})

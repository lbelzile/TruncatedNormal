test_that("trandn", {
  expect_error(trandn(0, -1));
  expect_error(trandn(c(0,1), c(1,2,3)));
  expect_vector(trandn(c(0,1,2,3), c(1,2,3,4)));
  expect_length(trandn(c(0,1), c(1,2)), 2);
  expect_lte(trandn(0,1), 1);
  expect_gte(trandn(0,1), 0);
})

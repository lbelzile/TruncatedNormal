test_that("ntail2", {
  expect_error(ntail2(c(1,1), c(2,2,2)));
  expect_error(ntail2(c(1, 3), c(2, 1)));
  expect_error(ntail2(c(0,0), c(2, 1)));
  expect_vector(ntail2(c(1,0.5), c(2,1)));
  expect_lte(ntail2(0.5,4), 4);
  expect_gte(ntail2(0.5,4), 0.5);
})

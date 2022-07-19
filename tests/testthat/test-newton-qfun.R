#source("../newton_qfun.cpp", chdir = TRUE)
library(testthat);

#library(TruncatedNormal);

#newton, qfun tester
test_that("newton3, qfun3",{
  expect_equal(newton3(c(0.5, 0.5), c(-1,-2), c(0, 1)), newton(c(0.5, 0.5), c(-1,-2), c(0, 1)));
  expect_length(newton3(c(0.5, 0.5), c(-1,-2), c(0, 1)), 2);
  expect_equal(qfun3(0.5), qfun(0.5));
  expect_equal(qfun3(-Inf), Inf);
  expect_equal(qfun3(Inf), 0);
  expect_length(qfun3(c(0,1,2,3,4)), 5);
  })


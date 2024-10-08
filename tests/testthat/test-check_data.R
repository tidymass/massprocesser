library(testthat)
test_that("check_data", {
  testthat::expect_error(check_targeted_table())
})

local_edition(3)

test_that("age_mat works", {
  expect_snapshot_value(age_mat(NS_params), style = "deparse")
})
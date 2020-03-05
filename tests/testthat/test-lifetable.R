context("lifetable")

test_that("Life-table is calculated correctly", {
  data(ovarian, package="survival")
  lt_df <- lifetable(ovarian, c(0, 59, 500, 1500), futime, fustat)
  expect_data_frame(lt_df, nrows = 3, ncols = 9)
  expect_equal(lt_df$n, c(26, 26, 12))
  expect_equal(lt_df$riskset, c(26, 24, 7))
  expect_equal(round(lt_df$hazard, 3), c(0, 0.417, 0.286))
  expect_equal(round(lt_df$survival, 3), c(1, 0.583, 0.417))

  # expect_data_frame(lifetable(ovarian, NULL, futime, fustat))

})

context("lwf_pec")

test_that("lwf_pec and pec give same results", {
  library(checkmate)

  library(eha)
  mod <- aftreg(Surv(enter, exit, event) ~ ses, data = mort)
  surv <- predictSurvProb(mod, newdata = mort[1:2, ], times = c(0, 5, 10))
  expect_matrix(surv, mode = "numeric", any.missing = FALSE,
    nrows = 2, ncols = 3)
  expect_identical(round(surv[1, ], 3), c(1.000, 0.962, 0.904))
  expect_identical(round(surv[2, ], 3), c(1.000, 0.940, 0.849))

})

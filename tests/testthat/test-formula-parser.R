context("formula-parser")

test_that("Formula creation works", {

  num_vars <- c("x")
  cat_vars <- c("z")
  form1 <- create_formulas(num_vars, cat_vars, lhs = "Surv(time, event)~")
  form2 <- create_formulas(num_vars, lhs = "Surv(time, event)~")
  form3 <- create_formulas(cat_vars, lhs = "Surv(time, event)~")
  form4 <- create_formulas(num_vars, cat_vars, lhs = "Surv(time, event)~",
    fname = "pspline")
  expect_identical(length(form1), 3L)
  expect_identical(length(form4), 5L)
  expect_identical(form4[3], "Surv(time, event)~pspline(x)+z")

})

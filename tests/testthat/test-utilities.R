context("Utility functions")
fit <- coxph(Surv(time, status)~pspline(karno, df=4), data=veteran)
term.karno <- get_terms(fit, veteran, terms="karno")


test_that("Interval infos correct", {
	expect_data_frame(int_df <- int_info(1:2), nrows = 2, ncols = 5)
	expect_equal(names(int_df), c("tstart", "tend", "intlen", "intmid", "interval"))
  expect_equal(int_df$interval, c("(0,1]", "(1,2]"))
	expect_equal(int_info(1:2, right_closed = FALSE)$interval, c("[0,1)", "[1,2)"))
})


test_that("Termplot extraction works correctly", {
  expect_equal(nrow(term.karno), 100)
  expect_equal(ncol(term.karno), 6)
  expect_equal(names(term.karno), c("term", "x", "eff", "se", "ci.lower", "ci.upper"))
})


test_that("Throws error on wrong input", {
  expect_error(get_terms(fit, veteran, terms = c("a")))
})

test_that("Split function works correctly", {
  data(ovarian, package="survival")
  split_df <- split_info(ovarian,"futime", "fustat",
    breaks = c(0, 59, 500, 1500))
  expect_data_frame(split_df, nrows = 63L, ncols = 6L)
  expect_subset(c("tstart", "tend", "interval", "status"), names(split_df))
  # expect_subset(c("cut", "id_var", "intvars"), names(attributes(split_df)))
  expect_equal(attr(split_df, "breaks"), c(0, 59, 500, 1500))

  split_df <- split_info(ovarian,"futime", "fustat",
    breaks = c(0, 59, 500, 1500), right_closed = FALSE)
  expect_data_frame(split_df, nrows = 64L, ncols = 6L)

})

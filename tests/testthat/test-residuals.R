context("Residuals")

## Cox-Snell Residuals
data("tongue", package="KMsurv")
cox.tongue <- coxph(Surv(time, delta)~type, data=tongue)

test_that("Dimensions and names correct", {
	expect_equal(length(get_csvec(cox.tongue)), nrow(tongue))
	expect_data_frame(get_coxsnell(cox.tongue), nrows=56L, ncols=3L)
	expect_equal(names(get_coxsnell(cox.tongue)),
		c("coxsnell", "cumu_hazard", "survival"))
})


test_that("Error thrown on misspecification", {
	expect_error(get_csvec())
	expect_error(get_csvec(iris))
	expect_error(get_coxsnell(iris))
})


test_that("Cox-Snell plots work", {
	expect_is(gg_coxsnell(cox.tongue, type="cumu_hazard"), c("gg", "ggplot"))
	expect_is(gg_coxsnell(cox.tongue, type="cdf"), c("gg", "ggplot"))
})


test_that("Cox-Snell throws error on misspecification", {
	expect_error(gg_coxsnell(iris, type="cumu_hazard"))
})


## Schoenfeld residuals
data("veteran", package="survival")

fit.vet <- coxph(Surv(time, status) ~ trt + celltype + karno+
	diagtime + age + prior, data=veteran)

zph.vet <- cox.zph(fit.vet, transform="identity", global=FALSE)

scaledsch <- get_scaledsch(fit.vet, transform="identity")


test_that("Scaled Schoenfeld residuals obtained correctly", {
	expect_equal(nrow(scaledsch), nrow(zph.vet$y)*ncol(zph.vet$y))
	expect_equal(dim(scaledsch), c(768, 4))
	expect_equal(names(scaledsch), c("time", "transform", "variable", "residual"))

})


test_that("Scaled Schoenfeld plots work", {
	expect_is(gg_scaledsch(fit.vet), c("gg", "ggplot"))
	expect_is(gg_scaledsch(fit.vet, transform="identity"), c("gg", "ggplot"))
})

test_that("Scaled Schoenfeld plot throws error on misspecification", {
	expect_error(gg_scaledsch(iris, type="cumu_hazard"))
})

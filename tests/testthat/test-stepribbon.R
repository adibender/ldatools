context("Step ribbon plots")

data("tongue", package="KMsurv")
library(survival)
cox.tongue <- coxph(Surv(time, delta)~type, data=tongue)
tdf <- broom::tidy(survfit(cox.tongue))

gg.obj <- ggplot(tdf, aes(x=time, y=estimate)) + 
	geom_stepribbon(aes(ymin=conf.low, ymax=conf.high))


test_that("Geom stepribbon works without error", {
	expect_is(gg.obj, c("gg", "ggplot"))
})
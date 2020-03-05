#' Overall model fit diagnostic via Cox-Snell residuals.
#'
#' \code{gg_coxsnell} extracts Cox-Snell residuals from an \code{\link[survival]{coxph}}
#' object (which for a correct model should be a censored sample from Exp(1)) and
#' the cumulative hazard rate of the residuals vs. the residuals.
#'
#' @inheritParams survival::cox.zph
#' @param type \code{character}. Optional argument. If \code{"cumu_hazard"} plots
#' cumulative hazard (Nelson-Aalen estimate) vs. Cox-Snell residuals.
#' If \code{"cdf"} plot empirical cumulative distribution Function (Breslow estimate)
#' vs. Cox-Snell residuals.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @examples
#' library(survival)
#' library(ggplot2)
#' data("tongue", package="KMsurv")
#' cox.tongue <- coxph(Surv(time, delta)~as.factor(type), data=tongue)
#' gg_coxsnell(cox.tongue) +
#'  geom_abline(intercept=0, slope=1, col=2)
#' gg_coxsnell(cox.tongue, type="cdf") +
#'  geom_line(aes(y=F), col=2)
#' @import dplyr ggplot2
#' @importFrom stats pexp
#' @importFrom rlang .data
#' @export
gg_coxsnell <- function(fit, type=c("cumu_hazard", "cdf")) {

  type <- match.arg(type)

  cs_data <- get_coxsnell(fit)

  if (type == "cdf") {
    cs_data <- cs_data %>%
      mutate(
        CDF = 1 - .data$survival,
        F   = pexp(.data$coxsnell))
    type      <- "CDF"
    ylab.text <- "Cumulative Distribution Function"
  } else {
    type      <- "cumu_hazard"
    ylab.text <- "Cumulative Hazard"
  }

  gg.coxsnell <- ggplot(cs_data, aes_string(x = "coxsnell", y = type)) +
    geom_point() +
    geom_step() +
    xlab("Cox-Snell residual") +
    ylab(ylab.text)

  return(gg.coxsnell)

}


#' Plot scaled schoenfeld residuals vs. time
#'
#' Creates a \code{ggplot} object/plot of scaled schoenfeld residuals
#' vs. (transformed) time.
#'
#' @inheritParams survival::cox.zph
#' @import ggplot2
#' @importFrom checkmate assert_class
#' @export
gg_scaledsch <- function(fit, transform = "km") {

  ## check input
  assert_class(fit, "coxph")

  ## obtain scaled schoenfeld residuals and (transformed time)
  scaledsch <- get_scaledsch(fit = fit, transform = transform)
  trans.string <- ifelse(unique(scaledsch$transform) == "identity", "t",
    paste0(unique(scaledsch$transform), "(t)"))

  gg.zph <- ggplot(scaledsch, aes_string(x = "time", y = "residual")) +
    geom_point() +
    facet_wrap(~variable, nrow = 2, scales = "free_y") +
    geom_smooth(method = "lm", lty = 2, aes(col = "lm")) +
    geom_smooth(method = "gam", formula = y~s(x), aes(col = "gam")) +
    scale_color_discrete(name = "method") +
    geom_hline(yintercept = 0, lty = 3) +
    xlab(trans.string) + ylab(expression(beta(t)))

  return(gg.zph)

}

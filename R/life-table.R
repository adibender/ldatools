#' Implementation of the actuarial life table method.
#'
#' @inherit split_info
#' @examples
#' data("Melanoma", package="MASS")
#' Melanoma$status <- ifelse(Melanoma$status == 1, 1, 0)
#' lifetable(Melanoma, breaks = seq(0, 6000, 1500))
#' @return A \code{\link[tibble]{tibble}} containing riskset information as well
#' as hazard and survival function calculated by using the actuarial life table
#' method.
#' @seealso split_info
#' @export
lifetable <- function(
  data,
  breaks       = NULL,
  time_var     = "time",
  status_var   = "status",
  right_closed = FALSE,
  max_end      = FALSE) {

  time_var   <- enquo(time_var)
  status_var <- enquo(status_var)

  if (is.null(breaks)) {
    breaks <- ceiling(seq(0, max(data[[quo_name(time_var)]]), length.out = 10))
  }

  split_df <- split_info(data = data, breaks = breaks,
    !!time_var, !!status_var, right_closed = right_closed,
    max_end = max_end) %>%
    select(-one_of("time")) %>%
    group_by(.data$id) %>%
    mutate(censored = 1 * (.data$status == 0 & max(.data$tend) == .data$tend)) %>%
    ungroup()

  split_df %>%
    group_by(.data$tstart, .data$tend, .data$ interval) %>%
    summarize(
      n        = n(),
      events   = sum(.data$status),
      dropouts = sum(.data$censored)) %>%
    ungroup() %>%
    mutate(
      riskset  = n - .data$dropouts / 2,
      hazard   = .data$events / .data$riskset,
      survival = cumprod(1 - .data$hazard)) %>%
    arrange(.data$tstart)

}

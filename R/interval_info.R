#' Create start/end times and interval information
#'
#' Given interval breaks points, returns data frame with information on
#' interval start time, interval end time, interval length and a interval factor
#' variable (left open intervals). If object of class ped is provided, extracts
#' unique interval information from object.
#'
#' @param x A numeric vector of cut points in which the follow-up should be
#' partitioned in or object of class \code{ped}.
#' @param ... Currently ignored.
#' @rdname int_info
#' @return A data frame containing the start and end times of the
#' intervals specified by the \code{x} argument. Additionally the interval
#' length, interval mid-point and a factor variable of the intervals themselves.
#' @export
int_info <- function(x, ...) {
  UseMethod("int_info",  x)
}


#' @inheritParams int_info
#' @inheritParams split_info
#' @param min_time Only intervals that have lower borders larger than
#' this value will be included in the resulting data frame.
#' @import checkmate dplyr tibble
#' @examples
#' ## create interval information from cut points
#' int_info(c(1, 2.3, 5))
#'
#' @rdname int_info
#' @export
int_info.default <- function(
  x,
  right_closed = TRUE,
  min_time   = 0L,
  ...) {

  # check inputs
  assert_numeric(x, lower = 0, any.missing = FALSE)
  assert_numeric(min_time, lower  = 0L)

  # sort x and add origin if necessary
  if (is.unsorted(x)) {
    x <- sort(x)
  }
  if (min(x != 0)) {
    x <- c(0, x)
  }

  intlen <- diff(x)
  tstart <- x[-length(x)]
  tend   <- tstart + intlen

  tdf <- tibble(
      tstart = tstart,
      tend   = tend,
      intlen = intlen) %>%
    mutate(
      intmid   = tstart + intlen / 2,
      interval = label_intervals(tstart, tend, right_closed = right_closed))

  filter(tdf, tstart >= min_time)

}


#' Information on intervals in which times fall
#'
#' @inheritParams int_info
#' @param x An object from which interval information can be obtained,
#' see \code{\link{int_info}}.
#' @param times A vector of times for which corresponding interval information
#' should be returned.
#' @param ... Further arguments passed to \code{\link[base]{findInterval}}.
#' @import dplyr
#' @return A \code{data.frame} containing information on intervals in which
#' values of \code{times} fall.
#' @examples
#' set.seed(111018)
#' brks <- c(0, 4.5, 5, 10, 30)
#' int_info(brks)
#' x <- runif(3, 0, 30)
#' x
#' get_intervals(brks, x)
#'
#' @seealso \code{\link[base]{findInterval}} \code{\link{int_info}}
#' @rdname get_intervals
#' @export
get_intervals <- function(x, times, ...) {
  UseMethod("get_intervals", x)
}

#' @inherit get_intervals
#' @inheritParams base::findInterval
#' @seealso \code{\link[base]{findInterval}}
#' @rdname get_intervals
#' @export
get_intervals.default <- function(
  x,
  times,
  left.open        = TRUE,
  rightmost.closed = TRUE,
  ...) {

  # check inputs
  assert_numeric(times, lower = 0, finite = TRUE, all.missing = FALSE)

  int_df <- int_info(x)
  int    <- findInterval(
    x                = times,
    vec              = union(int_df$tstart, int_df$tend),
    left.open        = left.open,
    rightmost.closed = rightmost.closed)

  int_df %>%
    slice(int) %>%
    mutate(times = times) %>%
    arrange(times) %>%
    select(times, everything())

}

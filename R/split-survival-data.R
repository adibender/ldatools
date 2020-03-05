#' Create index of breaks survived
#'
#' @inheritParams split_info
#' @param time The time of the event.
#' @keywords internal
survived_breaks <- function(
  time,
  breaks,
  right_closed   = TRUE) {

  if (!right_closed) {
    n_survived <- sum(time >= breaks)
  } else {
    n_survived <- sum(time > breaks)
  }

  seq_len(n_survived)

}

#' @keywords internal
interval_status <- function(
  time,
  tend,
  status,
  right_closed = TRUE) {

  if(!right_closed) {
    event <- time < tend & status == 1
  } else {
    event <- time <= tend & status == 1
  }

  event * 1

}


#' Create character vector of interval labels
#'
#' @inheritParams split_info
#' @param tstart Vector of interval start times.
#' @param tend Vector of interval end times.
#' @keywords internal
label_intervals <- function(tstart, tend, right_closed = TRUE) {

  if (right_closed) {
    left_bracket  <- "("
    right_closed_bracket <- "]"
  } else {
    left_bracket  <- "["
    right_closed_bracket <- ")"
  }
  paste0(left_bracket, tstart, ",", tend, right_closed_bracket)

}

#' @inheritParams split_info
#' @keywords internal
process_breaks <- function(data, breaks, time_var, status_var, max_end) {

  if (is.null(breaks)) {
    breaks <- unique(c(0, data[[time_var]][data[[status_var]] == 1]))
  }
  max_time <- max(max(data[[time_var]]), max(breaks))
  # sort interval break points in case they are not (so that interval factor
  # variables will be in correct ordering)
  breaks <- sort(breaks)

  # add last observation to breaks if necessary
  if (max_end & (max_time > max(breaks))) {
    breaks <- c(breaks, max_time)
  }

  breaks

}

#' @inheritParams split_info
#' @keywords internal
add_id <- function(data, id_var) {

    if(id_var %in% names(data)) {
      if (length(unique(data[[id_var]])) != nrow(data)) {
        stop(paste0("Specified ID variable (", id_var, ") must have same number of unique values as number of rows in 'data'."))
      }
    } else {
      data[[id_var]] <- seq_len(nrow(data))
    }

    data

}

#' Transform standard time-to-event data to interval data format.
#'
#' Given a data set in standard format with one row per observation unit,
#' transform the data set into interval data, where intervals are specified
#' via breaks. Each observation unit will have as many rows as interval
#' start points that it survived.
#'
#' @param data Data set from which time and status variables will be extracted.
#' @param time_var The variable storing event times.
#' @param status_var The variable storing the event indicator.
#' @param id_var The variable storing ID variable. Will be created if not
#' in data set.
#' @param breaks The time points of the interval borders.
#' @param right_closed Logical. If \code{TRUE} (default), intervals are assumed right_closed
#' closed and left open. If \code{FALSE} left closed and right_closed open.
#' @param max_end logical. Should the last interval span until the last
#' observed censoring or event time (if larger than the largest specified
#' cut point).
#' @import checkmate dplyr tibble
#' @importFrom purrr map map_dbl flatten_dbl
#' @examples
#' data(ovarian, package="survival")
#' split_info(ovarian, futime, fustat, breaks = c(0, 59, 500, 1500))
#' @export
split_info <- function(
  data,
  time_var     = "time",
  status_var   = "status",
  id_var       = "id",
  breaks       = NULL,
  right_closed = TRUE,
  max_end      = FALSE) {


  time_var   <- enquo(time_var)
  status_var <- enquo(status_var)
  id_var     <- enquo(id_var)

  assert_data_frame(data, min.rows = 2, min.cols = 2)
  assert_numeric(breaks, lower = 0, finite = TRUE, any.missing = FALSE,
    min.len = 1, null.ok = TRUE)
  assert_flag(right_closed)
  assert_flag(max_end)
  assert_subset(c(quo_name(time_var), quo_name(status_var)), names(data))

  breaks <- process_breaks(data, breaks, quo_name(time_var),
    quo_name(status_var), max_end)

  interval_df <- data %>%
    add_id(quo_name(id_var)) %>%
    select(one_of(
      c(quo_name(id_var), quo_name(time_var), quo_name(status_var)))) %>%
    rename(
      "time"   = !!time_var,
      "status" = !!status_var) %>%
    mutate(
      time = ifelse(.data$time > max(breaks), max(breaks), .data$time),
      status = ifelse(.data$time > max(breaks), 0, .data$status))
  n_survived <- map(data[[quo_name(time_var)]],
    ~survived_breaks(., breaks = breaks, right_closed = right_closed))
  last       <- map_dbl(n_survived, max)
  n_survived <- flatten_dbl(n_survived)
  slice_ind  <- rep(seq_len(nrow(data)), times = last)

  # create data frame with interval information
  int_df <- int_info(breaks, right_closed = right_closed) %>%
    select(one_of(c("tstart", "tend", "interval"))) %>%
    slice(n_survived)

  # final interval data that will be output
  interval_df <- interval_df %>%
    slice(slice_ind) %>%
    bind_cols(int_df) %>%
    mutate(status = interval_status(.data$time, .data$tend, .data$status,
      right_closed = right_closed))

  # set some attributes that might be usefull in other functions
  attr(interval_df, "id_var") <- quo_name(id_var)
  attr(interval_df, "breaks") <- breaks
  attr(interval_df, "right_closed")  <- right_closed
  attr(interval_df, "interval_vars") <-  c(quo_name(id_var), "tstart", "tend",
    "interval", "time", "status")

  # rearrange columns
  interval_df %>%
    select(one_of(attr(interval_df, "interval_vars")))

}

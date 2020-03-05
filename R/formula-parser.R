#' Wrapes variable in a function
#'
#' @param x A String the will be wrapped in \code{fname}.
#' @param fname A string of the special function in which \code{x} will be
#' wrapped in. Defaults to \code{"rcs"}
#' @seealso \code{\link[rms]{rms.trans}}
#' @keywords internal
fpaste <- function(x, fname="rcs") {
  paste0(fname, "(", x, ")")
}

#' Check if string appears in character vector multiple times
#'
#' @param char The string for which to check multiple appearances in \code{vec}.
#' @param vec The character vector in which to check multiple appearences of \code{char}.
#' @return \code{TRUE} if \code{char} appears more than once in \code{vec},
#' \code{FALSE} otherwise.
get_multappear <- function(char, vec) {
  sum(grepl(char, vec, fixed = TRUE)) > 1
}

#' Check multiple strings for multiple appearances
#'
#' A vectorized version of \code{get_multappear}.
#'
#' @param charvec A Vector of strings. Each will be checked for multiple appearences
#' in \code{vec}
#' @inheritParams get_multappear
#' @return \code{TRUE} if any element of \code{charvec} appears in \code{vec}
#' more than once, \code{FALSE} otherwise.
vget_multappear <- function(charvec, vec) {
  any(sapply(charvec, function(x) {
    get_multappear(x, vec)
  }))
}

#' Create formulas for all possible combinations of variables
#'
#' Creates formulas for all possible combinations of variables provided to
#' the function. All terms can appear as linear terms or non-linear terms.
#' @param num_vars A character vector of covariates on numeric scale.
#' @param cat_vars A character vector of categorical covariates
#' @param lhs The left-hand side of the formula to be produced.
#' @inheritParams fpaste
#' @importFrom e1071 bincombinations
#' @return A vector of formula strings for possible (and correct) model
#' specifications given the variables provided.
#' @examples
#' create_formulas(num_vars = c("x1", "x2"), lhs = "Surv(time, status) ~")
#' create_formulas(num_vars = c("x1", "x2"), cat_vars = c("z1"),
#'  lhs = "Surv(time, status) ~")
#' create_formulas(num_vars = c("x1"), cat_vars = c("z1", "pspline(x2)"),
#'  lhs = "Surv(time, status) ~")
#' create_formulas(num_vars = c("x1", "x2"), cat_vars = c("z1"),
#'  lhs = "Surv(time, status) ~", fname = "psline")
#' @export
create_formulas <- function(
  num_vars,
  cat_vars,
  lhs,
  fname = NULL) {

  if (missing(num_vars) & missing(cat_vars)) {
    stop("No variables specified")
  }

  if (!is.null(fname) & !is.null(num_vars)) {
    vars  <- c(num_vars, fpaste(num_vars, fname = fname))
  } else {
    vars <- num_vars
  }

  if (!missing(cat_vars)) {
    vars <- c(vars, cat_vars)
  }
  combos   <- bincombinations(length(vars))[-1, , drop = FALSE] == 1
  ind_mult <- apply(combos, 1, function(z)
    vget_multappear(num_vars, vars[z]))
  formxy   <- sapply(which(!ind_mult),
    function(z) {
        paste0(vars[combos[z, , drop = TRUE]], collapse = "+")
  })

  return(paste0(lhs, formxy))

}

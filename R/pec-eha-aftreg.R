#' Predict survival probabilities
#'
#' See \code{pec::predictSurvProb} for details.
#'
#' @name predictSurvProb
#' @rdname predictSurvProb
#' @keywords internal
#' @importFrom pec predictSurvProb
#' @export
NULL

#' Survival probability for eha::aftreg objects
#'
#' The function is named to be compatible with \code{pec} package. However,
#' \code{pec::pec} currently does not allow \code{method = 'aftreg'}. Therefore, this
#' function is only used internally in \code{pec_cv}.
#' The calculation of the Survival probability is based on the formula given
#' in \code{?eha::aftreg}.
#'
#' @inheritParams pec::predictSurvProb
#' @importFrom pec predictSurvProb
#' @importFrom stats model.matrix coef
#' @export
predictSurvProb.aftreg <- function(object, newdata, times, ...) {

  coefs_aft <- coef(object)
  ind_scale <- grep("log(scale)", names(coefs_aft), fixed = TRUE)
  ind_shape <- grep("log(shape)", names(coefs_aft), fixed = TRUE)

  coef_scale <- coefs_aft[ind_scale]
  coef_shape <- coefs_aft[ind_shape]
  beta       <- coefs_aft[c(-ind_scale, -ind_shape)]
  X    <- model.matrix(object, newdata)[, -1, drop = FALSE]
  etas <- as.numeric(X %*% beta)

  vapply(
    times,
    function(time) {
      exp(- (time / exp(coef_scale - etas) ) ** exp(coef_shape))
    },
    numeric(length(etas)))

}

#' Predict survival probabilities
#'
#' See \code{pec::ipcw} for details.
#'
#' @name ipcw
#' @rdname ipcw
#' @keywords internal
#' @importFrom pec ipcw
#' @export
NULL


#' Inverse probability censoring weights (IPCW) for left-truncated data
#'
#' This function uses \code{eha::aftreg} to estimate IPCW for left-truncated
#' data based on an AFT Weibull model. Initially designed to be used in
#' conjunction with \code{pec::pec}, however, \code{pec} does not support
#' \code{eha::aftreg} as a method for calculation of IPCWs. The function
#' is loosely based on \code{pec:::ipcw.cox} and only supports calculating
#' IPCWs at specified times (\code{what = "IPCW.times"}; see \code{what} argument
#' in \code{?pec::ipcw}).
#'
#' @inheritParams pec::ipcw
#' @importFrom pec ipcw
#' @importFrom stats na.omit
#' @export
ipcw.aftreg <- function(
  formula,
  data,
  method = "aftreg",
  args,
  times,
  subjectTimes,
  subjectTimesLag,
  what = "IPCW.times") {

    call <- match.call()
    #
    wdata               <- data
    status_var          <- all.vars(formula)[3]
    wdata[[status_var]] <- 1 - wdata[[status_var]]
    wform               <- formula
    stopifnot(NROW(na.omit(wdata)) > 0)

    if (missing(args) || is.null(args)) {
      args <- list(x = TRUE)
    }
    wdata$status[wdata$stopp > max(times)] <- 0
    wdata$stopp[wdata$stopp > max(times)] <- max(times)
    fit <- do.call(eha::aftreg, c(list(wform, data = wdata), args))
    #  weigths at requested times
    if (match("IPCW.times", what, nomatch = FALSE)) {
        IPCW.times <- predictSurvProb(fit, newdata = wdata, times = times,
          se.fit = FALSE)
    }
    else
        IPCW.times <- NULL
    #  weigths at subject specific event times
    if (match("IPCW.subjectTimes", what, nomatch = FALSE)){
        if (subjectTimesLag == 1)
            IPCW.subjectTimes <- predictSurvProb(fit, newdata = wdata,
              times = subjectTimes - min(diff(c(0, unique(subjectTimes)))) / 2)
        else if (subjectTimesLag == 0) {
            IPCW.subjectTimes <- predictSurvProb(fit, newdata = wdata,
              times = subjectTimes)
        }
        else stop("SubjectTimesLag must be 0 or 1")
    }
    else IPCW.subjectTimes <- NULL
    out <- list(
      times             = times,
      IPCW.times        = IPCW.times,
      IPCW.subjectTimes = IPCW.subjectTimes,
      fit               = fit,
      call              = call,
      method            = method)
    class(out) <- "IPCW"
    out
}

#' @keywords internal
has_censoring <- function(data, formula, times) {

  n_cens <- data %>%
    filter(.data$status == 0) %>%
    summarize(n_cens = sum(.data$stopp < max(.data$times))) %>%
    pull(.data$n_cens)

  n_cens > 0

}

#' Calculate prediction error curve
#' @import pec
#' @inheritParams pec::predictSurvProb
#' @inheritParams pec::ipcw
#' @param ipcw.form The formula used to calculate inverse probability censoring
#' weights.
#' @param ... Further arguments passed to \code{\link[pec]{predictSurvProb}}
#' @export
#' @examples
#' library(eha)
#' data(mort)
#' mort_train <- mort[1:1000,]
#' mort_test <- mort[-1:-1000,]
#' mod <- aftreg(Surv(enter, exit, event) ~ ses, param = "lifeExp",
#'  data = mort_train)
#' pec <- aft_pec(
#'  mod,
#'  ipcw.form = Surv(time = enter, time2=exit, event = event) ~ ses,
#'  times = seq(0, max(mort_test$exit) - 1, by = 1),
#'  newdata = mort_test)
#' plot(seq(0, max(mort_test$exit)-1, by = 1), pec, type = "l")
aft_pec <- function(
  object,
  ipcw.form = Surv(time = start, time2 = end, event = event) ~ 1,
  newdata,
  times,
  ...) {

  pec_out <- pec2(
    object     = object,
    formula    = ipcw.form,
    data       = newdata,
    times      = times,
    cens.model = "cox",
    exact      = FALSE,
    reference  = FALSE,
    ... )

  brier <- pec_out$AppErr[[1]]
  names(brier) <- paste0("t", times)
  attr(brier, "ibs") <- ibs(pec_out)

  return(brier)

}


#' Calculate cross-validated prediction error curves
#' @inheritParams pec::ipcw
#' @inheritParams modelr::crossv_kfold
#' @param data The data set used for cross-validated brier-score calculation.
#' @param mod.form The model specification used to fit the model in each fold.
#' Will be passed to \code{\link[eha]{aftreg}}.
#' @param ipcw.form The specification of the model that calculates the
#' inverse probability censoring weights (ipcw). Will be passed to
#' \code{\link[pec]{ipcw}}.
#' @param seed If specified this seed will be used (to replicate CV folds for example).
#' @param ... Further arguments passed to \code{\link{aft_pec}}.
#' @import survival
#' @importFrom modelr crossv_kfold
#' @importFrom stats as.formula
#' @importFrom eha aftreg
#' @importFrom dplyr as_tibble
#' @export
#' @examples
#' library(eha)
#' pec_cv_df <- aft_pec_cv(
#'  data = mort,
#'  mod.form = "Surv(time = enter, time2 = exit, event = event)~ses",
#'  ipcw.form = "Surv(time = enter, time2 = exit, event = event)~ses",
#'  times = seq(0, 10, by = 1),
#'  method = "cox",
#'  k = 5)
#' pec_cv_df <- aft_pec_cv(
#'  data = mort,
#'  mod.form = "Surv(time = enter, time2 = exit, event = event)~ses",
#'  ipcw.form = "Surv(time = enter, time2 = exit, event = event)~ses",
#'  times = seq(0, 10, by = 1),
#'  method = "aftreg",
#'  k = 5)
aft_pec_cv <- function(
  data,
  mod.form,
  ipcw.form,
  times  = c(0, 30, 100),
  seed   = NULL,
  k      = 10,
  ...) {

  if (!is.null(seed)) set.seed(seed)
  data_cv <- crossv_kfold(data, k)

  cv_folds <- list()
  for (i in seq_len(nrow(data_cv))) {
    # print(paste("Cross validation:", i))
    train  <- as_tibble(data_cv$train[[i]])
    test   <- as_tibble(data_cv$test[[i]])
    object <- aftreg(as.formula(mod.form), data = train)
    cv_folds[[i]] <- aft_pec(
      object    = object,
      ipcw.form = as.formula(ipcw.form),
      newdata   = test,
      times     = times,
      ...)
  }

  colMeans(do.call(rbind, cv_folds))

}

#' Wrapper functions that calculates CV prediction error
#'
#' @inheritParams aft_pec_cv
#' @inheritParams create_formulas
#' @param reference Logical. If \code{TRUE} (default) a baseline model
#' (\code{ ~ 1}) is added for comparison.
#' @importFrom tidyr gather nest
#' @import pbapply dplyr survival
#' @export
#' @examples
#' library(eha)
#' pec <- aft_pec_wrapper(
#'  data = mort,
#'  num_vars = "ses",
#'  ipcw.form = "Surv(time = enter, time2 = exit, event = event)~ses",
#'  times = seq(0, 17, by = 1),
#'  k = 5,
#'  lhs = "Surv(time = enter, time2 = exit, event = event)~",
#'  fname = "")
aft_pec_wrapper <- function(
  data,
  num_vars,
  cat_vars,
  ipcw.form,
  times,
  lhs      = "Surv(time = start, time2 = end, event = event)~",
  seed     = NULL,
  k        = 10,
  fname    = NULL,
  reference = TRUE,
  ...) {

  forms <- create_formulas(num_vars, cat_vars, lhs = lhs, fname = fname)
  if(reference) {
    forms <- c(paste0(lhs, "1"), forms)
  }
  pec.all <- pblapply(forms,
    function(z) {
      aft_pec_cv(
        data      = data,
        mod.form  = as.formula(z),
        ipcw.form = ipcw.form,
        times     = times,
        k         = k)
  })

  pec_df <- do.call(rbind, pec.all) %>% as_tibble() %>%
    mutate(
      id      = factor(row_number()),
      formula = forms) %>%
    gather("time", "brier", one_of(paste0("t", times))) %>%
    mutate(time2 = as.numeric(sub("t", "", .data$time))) %>%
    nest(-.data$id, -.data$formula, .key = "pec") %>%
    select(id, pec, everything())

  return(pec_df)

}


#' Calculate integrated brier score
#'
#' @param pec_df A data frame containing prediction error (brier score)
#' on different time points.
#' @return A data frame containing the integrated brier socre (IBS) for
#' each id. Results are sorted with respect to lowest IBS.
#' @importFrom tidyr unnest
#' @export
integrated_brier <- function(pec_df) {

  pec_df %>%
    unnest() %>%
    filter(!is.na(.data$brier)) %>%
    group_by(.data$id, .data$formula) %>%
    # see pec:::cprs, especially pec:::Dint
    summarize(ibs = sum(.data$brier * c(0, diff(.data$time2))) /
      (max(.data$time2) - min(.data$time2))) %>%
    arrange(.data$ibs) %>%
    ungroup()

}

#' Extract formula of the best model
#'
#' @inheritParams integrated_brier
#' @param times The time(s) for which the best fit will be searched.
#' @import dplyr
#' @export
#'
get_best <- function(
  pec_df,
  times = NULL) {

  if (is.null(times)) {
    integrated_brier(pec_df) %>%
      filter(.data$ibs == min(.data$ibs))
  } else {
    pec_df %>%
      unnest() %>%
      filter(.data$time2 %in% times) %>%
      group_by(.data$time2) %>%
      filter(.data$brier == min(.data$brier)) %>%
      ungroup()
  }
}

#' Extract specification of model with the lowest brier score
#'
#' @inheritParams get_best
#' @return The formula of the model that produced the lowest brier score at
#' time \code{time}.
#' @importFrom tidyr unnest
#' @export
get_best_form <- function(pec_df, times = NULL) {

  if (is.null(times)) {
    pec_df %>%
      integrated_brier() %>%
      filter(.data$ibs == min(.data$ibs)) %>%
      select(one_of("formula")) %>% unnest()
  } else {
    pec_df %>%
      unnest() %>%
      filter(.data$time2 %in% times) %>%
      filter(.data$brier == min(.data$brier)) %>%
      select(one_of("formula")) %>% unnest()
  }

}

#' Clone of pec::pec function
#'
#' Slightly modified such that start/stop formulas are accepted and processed.
#' @inheritParams pec::pec
#' @importFrom stats model.frame model.response update terms median na.fail
#' @keywords internal
pec2 <- function (object, formula, data, traindata, times, cause, start,
    maxtime, exact = TRUE, exactness = 100, fillChar = NA, cens.model = "cox",
    ipcw.refit = FALSE, ipcw.args = NULL, splitMethod = "none",
    B, M, reference = TRUE, model.args = NULL, model.parms = NULL,
    keep.index = FALSE, keep.matrix = FALSE, keep.models = FALSE,
    keep.residuals = FALSE, keep.pvalues = FALSE, noinf.permute = FALSE,
    multiSplitTest = FALSE, testIBS, testTimes, confInt = FALSE,
    confLevel = 0.95, verbose = TRUE, savePath = NULL, slaveseed = NULL,
    na.action = na.fail, ...)
{
    theCall = match.call()
    if (match("replan", names(theCall), nomatch = FALSE))
        stop("The argument name 'replan' has been replaced by 'splitMethod'.")
    if (!missing(testIBS) && (!(is.logical(testIBS) || (length(testIBS) ==
        2 && is.numeric(testIBS)))))
        stop("Argument testIBS can be TRUE/FALSE or a vector of two numeric values.")
    if (missing(testIBS))
        testIBS <- FALSE
    if (keep.residuals && missing(testTimes))
        stop("To keep.residuals please specify testTimes.")
    if (missing(splitMethod) && multiSplitTest == TRUE) {
        stop("Need data splitting to compute van de Wiel's test")
    }
    if (missing(M) && multiSplitTest)
        M <- NA
    if (class(object)[1] != "list") {
        object <- list(object)
    }
    if (missing(formula)) {
        if (length(grep("~", as.character(object[[1]]$call$formula))) ==
            0) {
            stop(paste("Argument formula is missing and first model has no usable formula:",
                as.character(object[[1]]$call$formula)))
        }
        else {
            ftry <- try(formula <- eval(object[[1]]$call$formula),
                silent = TRUE)
            if ((class(ftry) == "try-error") || match("formula",
                class(formula), nomatch = 0) == 0)
                stop("Argument formula is missing and first model has no usable formula.")
            else if (verbose)
                warning("Formula missing. Using formula from first model")
        }
    }
    formula.names <- try(all.names(formula), silent = TRUE)
    if (!(formula.names[1] == "~") || (match("$", formula.names,
        nomatch = 0) + match("[", formula.names, nomatch = 0) >
        0)) {
        stop("Invalid specification of formula.\n Could be that you forgot the right hand side:\n ~covariate1 + covariate2 + ...?\nNote that any subsetting, ie data$var or data[,\"var\"], is not supported by this function.")
    }
    else {
        if (!(formula.names[2] %in% c("Surv", "Hist")))
            survp <- FALSE
        else survp <- TRUE
    }
    if (missing(data)) {
        data <- eval(object[[1]]$call$data)
        if (match("data.frame", class(data), nomatch = 0) ==
            0)
            stop("Argument data is missing.")
        else if (verbose)
            warning("Argument data is missing. I use the data from the call to the first model instead.")
    }
    cens.model <- match.arg(cens.model, c("cox", "marginal",
        "nonpar", "aalen", "none", "rfsrc", "aftreg"))
    histformula <- formula
    m <- model.frame(histformula, data, na.action = na.action)
    response <- model.response(m)
    if (match("Surv", class(response), nomatch = 0) != 0) {
        attr(response, "model") <- "survival"
        attr(response, "cens.type") <- "rightCensored"
        attr(response, "type") <- "counting"
        model.type <- "survival"
    }
    model.type <- attr(response, "model")
    predictHandlerFun <- "predictSurvProb"
    if (reference == TRUE) {
        ref_form <- as.formula(update(formula, ".~NULL"))
        ref_fit <- eha::aftreg(formula = ref_form, data = data, x = TRUE)
        ref_fit$call$data <- as.character(substitute(data))
        ref_fit$call$formula = ref_form
        ref_fit$formula <- as.formula(ref_fit$formula)
        object <- c(list(Reference = ref_fit), object)
    }
    if (is.null(names(object))) {
        names(object) <- sapply(object, function(o) class(o)[1])
        names(object) <- make.names(names(object), unique = TRUE)
    }
    else {
        if (any(names(object) == "")) {
            names(object)[(names(object) == "")] <- sapply(object[(names(object) ==
                "")], function(o) class(o)[1])
            names(object) <- make.names(names(object), unique = TRUE)
        }
        else {
        }
    }
    NF <- length(object)
    if (survp) {
        neworder <- order(response[, "stop"], -response[, "status"])
        if (predictHandlerFun == "predictEventProb") {
            event <- prodlim::getEvent(response, mode = "character")
            event <- event[neworder]
        }
        response <- response[neworder, , drop = FALSE]
        Y <- response[, "stop"]
        status <- response[, "status"]
    }
    else {
        cens.model <- "none"
        neworder <- order(response)
        Y <- response[neworder]
        status <- rep(1, length(Y))
    }
    if (predictHandlerFun == "predictEventProb") {
        availableCauses <- unique(event)
        if (!match(cause, availableCauses, nomatch = FALSE))
            stop("Cause ", cause, " is not among the available causes: ",
                paste(availableCauses, collapse = ", "))
        event <- event == cause
    }
    data <- data[neworder, ]
    unique.Y <- unique(Y)
    N <- length(Y)
    NU <- length(unique.Y)
    splitMethod <- resolvesplitMethod(splitMethod = splitMethod,
        B = B, N = N, M = M)
    B <- splitMethod$B
    ResampleIndex <- splitMethod$index
    k <- splitMethod$k
    do.resample <- !(is.null(ResampleIndex))
    if (keep.matrix == TRUE & !do.resample) {
        warning("Argument keep.matrix set to FALSE, since no resampling/crossvalidation is requested.")
        keep.matrix <- FALSE
    }
    if (missing(maxtime) || is.null(maxtime))
        maxtime <- unique.Y[NU]
    if (missing(start))
        if (survp == TRUE)
            start <- 0
        else start <- min(unique.Y)
    if (missing(times)) {
        if (exact == TRUE)
            times <- unique(c(start, unique.Y))
        else times <- seq(start, maxtime, (maxtime - start)/exactness)
    }
    else {
        if (exact == TRUE)
            times <- sort(c(start, unique(times), unique.Y))
        else times <- sort(unique(c(start, times)))
    }
    times <- times[times <= maxtime]
    NT <- length(times)
    if ((cens.model %in% c("aalen", "cox", "nonpar"))) {
        if (all(as.numeric(status) == 1) || sum(status) == N) {
            if (verbose)
                message("No censored observations: cens.model coerced to \"none\".")
            cens.model <- "none"
        }
        if ((cens.model != "nonpar") && length(attr(terms(formula),
            "factors")) == 0) {
            if (verbose == TRUE)
                message("No covariates  specified: Kaplan-Meier for censoring times used for weighting.")
            cens.model <- "marginal"
        }
    }
    if (predictHandlerFun == "predictEventProb") {
        iFormula <- as.formula(paste("Surv(itime,istatus)", "~",
            as.character(formula)[[3]]))
        iData <- data
        iData$itime <- response[, "stop"]
        iData$istatus <- response[, "status"]
        if (ipcw.refit == TRUE)
            stop("pec: internal refitting of censoring distribution not (not yet) supported for competing risks")
        ipcw.call <- NULL
        ipcw <- ipcw(formula = iFormula, data = iData, method = cens.model,
            args = ipcw.args, times = times, subjectTimes = Y,
            subjectTimesLag = 1, what = "IPCW.times")
        ipcw$dim <- if (cens.model %in% c("marginal", "none"))
            0
        else 1
    } else {
        if (ipcw.refit == TRUE && splitMethod$internal.name %in%
            c("Boot632plus", "BootCv", "Boot632"))
            ipcw.call <- list(formula = formula, data = NULL,
                method = cens.model, times = times, subjectTimes = NULL,
                subjectTimesLag = 1)
        else ipcw.call <- NULL
        ipcw <- ipcw2(formula = formula, data = data, method = cens.model,
            args = ipcw.args, times = times, subjectTimes = Y,
            subjectTimesLag = 1)
        ipcw$dim <- if (cens.model %in% c("marginal", "none"))
            0
        else 1
    }
    if (do.resample) {
        cm <- pec:::checkModels(object = object, model.args = model.args,
            model.parms = model.parms, splitMethod = splitMethod$internal.name)
        model.args <- cm$model.args
        model.parms <- cm$model.parms
    }
    AppErr <- lapply(1:NF, function(f) {
        fit <- object[[f]]
        extraArgs <- model.args[[f]]
        if (predictHandlerFun == "predictEventProb") {
            pred <- do.call(predictHandlerFun, c(list(object = fit,
                newdata = data, times = times, cause = cause),
                extraArgs))
            if (class(fit)[[1]] %in% c("matrix", "numeric"))
                pred <- pred[neworder, , drop = FALSE]
            .C("pecCR", pec = double(NT), as.double(Y), as.double(status),
                as.double(event), as.double(times), as.double(pred),
                as.double(ipcw$IPCW.times), as.double(ipcw$IPCW.subjectTimes),
                as.integer(N), as.integer(NT), as.integer(ipcw$dim),
                as.integer(is.null(dim(pred))), NAOK = TRUE,
                PACKAGE = "pec")$pec
        }
        else {
            pred <- do.call(predictHandlerFun, c(list(object = fit,
                newdata = data, times = times), extraArgs))
            if (class(fit)[[1]] %in% c("matrix", "numeric"))
                pred <- pred[neworder, , drop = FALSE]
            .C("pecSRC", pec = double(NT), as.double(Y), as.double(status),
                as.double(times), as.double(pred), as.double(ipcw$IPCW.times),
                as.double(ipcw$IPCW.subjectTimes), as.integer(N),
                as.integer(NT), as.integer(ipcw$dim), as.integer(is.null(dim(pred))),
                NAOK = TRUE, PACKAGE = "pec")$pec
        }
    })
    names(AppErr) <- names(object)
    if (splitMethod$internal.name %in% c("Boot632plus")) {
        if (verbose == TRUE) {
            message("Computing noinformation error using all permutations")
        }
        if (noinf.permute == FALSE) {
            NoInfErr <- lapply(1:NF, function(f) {
                fit <- object[[f]]
                extraArgs <- model.args[[f]]
                if (predictHandlerFun == "predictEventProb") {
                  pred <- do.call(predictHandlerFun, c(list(object = fit,
                    newdata = data, times = times, cause = cause),
                    extraArgs))
                }
                else {
                  pred <- do.call(predictHandlerFun, c(list(object = fit,
                    newdata = data, times = times), extraArgs))
                }
                if (predictHandlerFun == "predictEventProb")
                  .C("pec_noinfCR", pec = double(NT), as.double(Y),
                    as.double(status), as.double(event), as.double(times),
                    as.double(pred), as.double(ipcw$IPCW.times),
                    as.double(ipcw$IPCW.subjectTimes), as.integer(N),
                    as.integer(NT), as.integer(ipcw$dim), as.integer(is.null(dim(pred))),
                    NAOK = TRUE, PACKAGE = "pec")$pec
                else .C("pec_noinf", pec = double(NT), as.double(Y),
                  as.double(status), as.double(times), as.double(pred),
                  as.double(ipcw$IPCW.times), as.double(ipcw$IPCW.subjectTimes),
                  as.integer(N), as.integer(NT), as.integer(ipcw$dim),
                  as.integer(is.null(dim(pred))), NAOK = TRUE,
                  PACKAGE = "pec")$pec
            })
            names(NoInfErr) <- names(object)
        }
        else {
            if (verbose == TRUE) {
                message("Noinformation error simulation loop (B=",
                  B, ")")
            }
            NoInfErrList <- lapply(1:B, function(b) {
                if (verbose == TRUE) {
                  pec:::internalTalk(b, B, sign = ".")
                }
                responseNames <- colnames(response)
                noinf.b <- data[sample(1:NROW(data), replace = FALSE),
                  -match(responseNames, names(data))]
                noinf.b[, responseNames] <- response
                ipcw.b <- ipcw(formula = formula, data = noinf.b,
                  method = cens.model, args = ipcw.args, times = times,
                  subjectTimes = Y, subjectTimesLag = 1)
                noinfPredErr <- lapply(1:NF, function(f) {
                  fit.b <- pec:::internalReevalFit(object = object[[f]],
                    data = noinf.b, step = b, silent = FALSE,
                    verbose = verbose)
                  extraArgs <- model.args[[f]]
                  pred.b <- do.call(predictHandlerFun, c(list(object = fit.b,
                    newdata = noinf.b, times = times), extraArgs))
                  if (predictHandlerFun == "predictEventProb") {
                    pred.b <- do.call(predictHandlerFun, c(list(object = fit.b,
                      newdata = noinf.b, times = times, cause = cause),
                      extraArgs))
                    .C("pecCR", pec = double(NT), as.double(Y),
                      as.double(status), as.double(event), as.double(times),
                      as.double(pred.b), as.double(ipcw.b$IPCW.times),
                      as.double(ipcw.b$IPCW.subjectTimes), as.integer(N),
                      as.integer(NT), as.integer(ipcw$dim), as.integer(is.null(dim(pred.b))),
                      NAOK = TRUE, PACKAGE = "pec")$pec
                  }
                  else {
                    pred.b <- do.call(predictHandlerFun, c(list(object = fit.b,
                      newdata = noinf.b, times = times), extraArgs))
                    .C("pecSRC", pec = double(NT), as.double(Y),
                      as.double(status), as.double(times), as.double(pred.b),
                      as.double(ipcw.b$IPCW.times), as.double(ipcw.b$IPCW.subjectTimes),
                      as.integer(N), as.integer(NT), as.integer(ipcw$dim),
                      as.integer(is.null(dim(pred.b))), NAOK = TRUE,
                      PACKAGE = "pec")$pec
                  }
                })
                noinfPredErr
            })
            NoInfErrMat <- lapply(1:NF, function(f) {
                do.call("rbind", lapply(NoInfErrList, function(x) {
                  x[[f]]
                }))
            })
            NoInfErr <- lapply(NoInfErrMat, colMeans)
            names(NoInfErr) <- names(object)
        }
    }
    if (splitMethod$internal.name %in% c("crossval", "loocv")) {
        kCV <- pec:::kFoldCrossValidation(object = object, data = data,
            Y = Y, status = status, event = event, times = times,
            cause = cause, ipcw = ipcw, splitMethod = splitMethod,
            giveToModel = model.args, predictHandlerFun = predictHandlerFun,
            keep = keep.matrix, verbose = verbose)
        CrossValErr <- kCV$CrossValErr
        if (keep.matrix && B > 1)
            CrossValErrMat <- kCV$CrossValErrMat
    }
    if (splitMethod$internal.name %in% c("Boot632plus", "BootCv",
        "Boot632")) {
        if (verbose == TRUE) {
            message("Split sample loop (B=", B, ")")
        }
        if (missing(testTimes)) {
            testTimes <- NULL
        }
        BootCv <- pec:::bootstrapCrossValidation(object = object, data = data,
            Y = Y, status = status, event = event, times = times,
            cause = cause, ipcw = ipcw, ipcw.refit = ipcw.refit,
            ipcw.call = ipcw.call, splitMethod = splitMethod,
            multiSplitTest = multiSplitTest, testIBS = testIBS,
            testTimes = testTimes, confInt = confInt, confLevel = confLevel,
            getFromModel = model.parms, giveToModel = model.args,
            predictHandlerFun = predictHandlerFun, keepMatrix = keep.matrix,
            keepResiduals = keep.residuals, verbose = verbose,
            savePath = savePath, slaveseed = slaveseed)
        BootstrapCrossValErr <- BootCv$BootstrapCrossValErr
        Residuals <- BootCv$Residuals
        names(BootstrapCrossValErr) <- names(object)
        if (multiSplitTest == TRUE) {
            comparisons <- pec:::allComparisons(names(object))
            multiSplitTestResults <- list(testIBS = testIBS,
                B = B, M = M, N = N, testTimes = testTimes)
            multiSplitTestResults$Comparisons <- lapply(1:length(comparisons),
                function(cc) {
                  if (length(testTimes) > 0) {
                    allPairwisePvaluesTimes <- do.call("rbind",
                      lapply(BootCv$testedResid, function(b) {
                        b$pValue[[cc]]
                      }))
                    out <- list(pValueTimes = apply(allPairwisePvaluesTimes,
                      2, median))
                    if (keep.pvalues == TRUE) {
                      out$allPairwisePvaluesTimes <- allPairwisePvaluesTimes
                    }
                  }
                  else out <- NULL
                  if (length(testIBS) > 0) {
                    allPairwisePvaluesIBS <- sapply(BootCv$testedResid,
                      function(b) {
                        b$IBSpValue[[cc]]
                      })
                    out$pValueIBS <- median(allPairwisePvaluesIBS)
                  }
                  if (keep.pvalues == TRUE) {
                    out$allPairwisePvaluesIBS <- allPairwisePvaluesIBS
                  }
                  out
                })
            names(multiSplitTestResults$Comparisons) <- names(comparisons)
            class(multiSplitTestResults) <- "multiSplitTest"
        }
        if (keep.matrix == TRUE) {
            BootstrapCrossValErrMat <- BootCv$BootstrapCrossValErrMat
            names(BootstrapCrossValErr) <- names(object)
        }
    }
    if (splitMethod$internal.name == "Boot632") {
        B632Err <- lapply(1:NF, function(f) {
            0.368 * AppErr[[f]] + 0.632 * BootstrapCrossValErr[[f]]
        })
        names(B632Err) <- names(object)
    }
    if (splitMethod$internal.name == "Boot632plus") {
        B632plusErr <- lapply(1:NF, function(f) {
            Err1 <- pmin(BootstrapCrossValErr[[f]], NoInfErr[[f]])
            overfit <- (Err1 - AppErr[[f]])/(NoInfErr[[f]] -
                AppErr[[f]])
            overfit[!(Err1 > AppErr[[f]])] <- 0
            w <- 0.632/(1 - 0.368 * overfit)
            B632plusErr <- (1 - w) * AppErr[[f]] + w * Err1
            B632plusErr
        })
        names(B632plusErr) <- names(object)
    }
    out <- switch(splitMethod$internal.name, noPlan = list(AppErr = AppErr),
        Boot632plus = list(AppErr = AppErr, BootCvErr = BootstrapCrossValErr,
            NoInfErr = NoInfErr, Boot632plusErr = B632plusErr),
        Boot632 = list(AppErr = AppErr, BootCvErr = BootstrapCrossValErr,
            Boot632Err = B632Err), BootCv = list(AppErr = AppErr,
            BootCvErr = BootstrapCrossValErr), loocv = list(AppErr = AppErr,
            loocvErr = CrossValErr), crossval = list(AppErr = AppErr,
            crossvalErr = CrossValErr), noinf = list(AppErr = AppErr,
            NoInfErr = NoInfErr))
    observed.maxtime <- sapply(out, function(x) {
        lapply(x, function(y) {
            times[length(y) - sum(is.na(y))]
        })
    })
    minmaxtime <- min(unlist(observed.maxtime))
    if (multiSplitTest == TRUE) {
        out <- c(out, list(multiSplitTest = multiSplitTestResults))
    }
    if (keep.residuals == TRUE) {
        out <- c(out, list(Residuals = Residuals))
    }
    if (keep.matrix == TRUE && splitMethod$internal.name != "noPlan") {
        if (splitMethod$internal.name %in% c("crossval", "loocv")) {
            if (B > 1)
                out <- c(out, list(CrossValErrMat = CrossValErrMat))
        }
        else {
            if (splitMethod$internal.name != "noinf")
                out <- c(out, list(BootstrapCrossValErrMat = BootstrapCrossValErrMat))
        }
    }
    if (!is.na(fillChar))
        out <- lapply(out, function(o) {
            o[is.na(o)] <- fillChar
            o
        })
    if (!is.null(model.parms))
        out <- c(out, list(ModelParameters = BootCv$ModelParameters))
    if (!keep.index)
        splitMethod$index <- NULL
    n.risk <- N - prodlim::sindex(Y, times)
    if (keep.models == TRUE) {
        outmodels <- object
    }
    else {
        outmodels <- names(object)
        names(outmodels) <- names(object)
    }
    out <- c(out, list(call = theCall, response = model.response(m),
        time = times, n.risk = n.risk, models = outmodels, maxtime = maxtime,
        observed.maxtime = observed.maxtime, minmaxtime = minmaxtime,
        reference = reference, start = min(times), cens.model = cens.model,
        exact = exact, splitMethod = splitMethod))
    class(out) <- "pec"
    out
}


#' @inherit pec::ipcw
ipcw2 <- function (formula, data, method, args, times, subjectTimes,
  subjectTimesLag = 1, what)
{
    if (!missing(what))
        stopifnot(all(match(what, c("IPCW.times", "IPCW.subjectTimes"))))
    if (missing(what) || match("IPCW.times", what, nomatch = FALSE)) {
        stopifnot(length(times) > 0)
    }
    if(method == "cox") {
      method <- "cox2"
    }
    class(method) <- method
    UseMethod("ipcw", method)
}

#' @inherit pec::ipcw.cox
#' @keywords internal
ipcw.cox2 <- function (formula, data, method, args, times, subjectTimes,
  subjectTimesLag, what)
{
    if (missing(subjectTimesLag))
        subjectTimesLag = 1
    if (missing(what))
        what = c("IPCW.times", "IPCW.subjectTimes")
    call <- match.call()
    EHF <- prodlim::EventHistory.frame(formula, data, specials = c("strat"),
        stripSpecials = c("strat"), specialsDesign = FALSE, unspecialsDesign = FALSE)
    if (is.null(EHF$strat))
        wdata <- data.frame(cbind(unclass(EHF$event.history),
            EHF$design))
    else wdata <- data.frame(cbind(unclass(EHF$event.history),
        EHF$design, EHF$strat))
    wdata$status <- 1 - wdata$status
    wform <- update(formula, "Surv(entry, time, status)~.")
    stopifnot(NROW(na.omit(wdata)) > 0)
    if (missing(args) || is.null(args))
        args <- list(x = TRUE, y = TRUE, eps = 1e-06)
    args$surv <- TRUE
    fit <- do.call(rms::cph, c(list(wform, data = wdata), args))
    if (match("IPCW.times", what, nomatch = FALSE)) {
        IPCW.times <- rms::survest(fit, newdata = wdata, times = times,
            se.fit = FALSE)$surv
    }
    else IPCW.times <- NULL
    if (match("IPCW.subjectTimes", what, nomatch = FALSE)) {
        if (subjectTimesLag == 1)
            IPCW.subjectTimes <- rms::survest(fit, times = subjectTimes -
                min(diff(c(0, unique(subjectTimes))))/2, what = "parallel")
        else if (subjectTimesLag == 0) {
            IPCW.subjectTimes <- rms::survest(fit, times = subjectTimes,
                what = "parallel")
        }
        else stop("SubjectTimesLag must be 0 or 1")
    }
    else IPCW.subjectTimes <- NULL
    out <- list(times = times, IPCW.times = IPCW.times, IPCW.subjectTimes = IPCW.subjectTimes,
        fit = fit, call = call, method = method)
    class(out) <- "IPCW"
    out
}

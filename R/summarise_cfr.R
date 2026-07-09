#' Naive deaths / cases ratio
#'
#' @param n_deaths,n_cases Death and case counts.
#' @return The naive ratio, or `NA` when there are no cases.
#' @noRd
naive_cfr <- function(n_deaths, n_cases) {
  if (n_cases == 0) NA_real_ else n_deaths / n_cases
}

#' Posterior draws of the CFR and the onset-to-death delay moments
#'
#' Transforms the raw `brms` draws into the reported quantities: `cfr` (from the
#' logit-scale `cfr` intercept) and the delay `mean`/`sd` in days, computed from
#' the family's native parameters. Handles the lognormal and gamma families.
#' @param object A `cfrnow_fit`.
#' @return A `draws_df` with columns `cfr`, `delay_mean`, `delay_sd`.
#' @noRd
.cfr_quantities <- function(object) {
  dr <- posterior::as_draws_df(object)
  if (!"b_cfr_Intercept" %in% posterior::variables(dr)) {
    stop("summary() supports intercept-only `cfr` fits; for `cfr ~ covariates` ",
         "use brms::posterior_epred() on the fit directly.", call. = FALSE)
  }
  cfr <- stats::plogis(dr[["b_cfr_Intercept"]])
  if (object$cfrnow$family == "lognormal") {
    meanlog <- dr[["b_Intercept"]]                     # identity link
    sdlog <- exp(dr[["b_sigma_Intercept"]])            # log link
    delay_mean <- exp(meanlog + sdlog^2 / 2)
    delay_sd <- sqrt(exp(sdlog^2) - 1) * delay_mean
  } else {                                             # gamma
    delay_mean <- exp(dr[["b_Intercept"]])             # log link on the mean
    shape <- exp(dr[["b_shape_Intercept"]])            # log link
    delay_sd <- delay_mean / sqrt(shape)
  }
  res <- data.frame(cfr = cfr, delay_mean = delay_mean, delay_sd = delay_sd,
                    .chain = dr$.chain, .iteration = dr$.iteration,
                    .draw = dr$.draw)
  posterior::as_draws_df(res)
}

#' Summarise a mixture-cure CFR fit
#'
#' Reports the corrected CFR and the onset-to-death delay (mean and sd, in days)
#' as posterior quantiles with convergence diagnostics (`rhat`, `ess_bulk`). The
#' naive `deaths / cases` ratio is returned as an attribute; in real time it
#' underestimates the corrected CFR because not every fatal case has died by the
#' cut-off.
#'
#' When few deaths have resolved (a young outbreak) the CFR is only weakly
#' identified and its posterior stays close to the prior. This is reported via
#' the `cfr_low_information` attribute: `TRUE` when the CFR posterior sd exceeds
#' `info_tol` times the prior sd.
#'
#' @param object A `cfrnow_fit` from [fit_cfr()].
#' @param probs Quantiles to report.
#' @param info_tol Low-information threshold: flag when the CFR posterior sd is
#'   more than this fraction of the prior sd. Defaults to 0.9.
#' @param ... Unused.
#' @return A data frame with one row per quantity (`cfr`, `delay_mean`,
#'   `delay_sd`), carrying `naive_cfr`, `n_cases`, `n_deaths`, `cfr_prior_sd`
#'   and `cfr_low_information` attributes.
#' @family fit
#' @export
summary.cfrnow_fit <- function(object, probs = c(0.025, 0.5, 0.975),
                               info_tol = 0.9, ...) {
  if (!inherits(object, "cfrnow_fit")) {
    stop("`object` must come from fit_cfr().", call. = FALSE)
  }
  q <- .cfr_quantities(object)
  qcols <- paste0("q", probs * 100)
  sm <- posterior::summarise_draws(
    q, mean = mean,
    stats::setNames(lapply(probs, function(p) {
      function(x) stats::quantile(x, p, names = FALSE)
    }), qcols),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  )
  out <- as.data.frame(sm)
  names(out)[1] <- "quantity"

  cfr_post_sd <- stats::sd(posterior::extract_variable(q, "cfr"))
  attr(out, "naive_cfr") <- naive_cfr(object$cfrnow$n_deaths,
                                      object$cfrnow$n_cases)
  attr(out, "n_cases") <- object$cfrnow$n_cases
  attr(out, "n_deaths") <- object$cfrnow$n_deaths
  attr(out, "cfr_prior_sd") <- object$cfrnow$cfr_prior_sd
  attr(out, "cfr_low_information") <-
    !is.na(object$cfrnow$cfr_prior_sd) &&
      cfr_post_sd > info_tol * object$cfrnow$cfr_prior_sd
  out
}

#' @rdname summary.cfrnow_fit
#' @param x A `cfrnow_fit`.
#' @export
print.cfrnow_fit <- function(x, ...) {
  s <- summary(x)
  message("<cfrnow_fit> ", x$cfrnow$family, " delay")
  message("  cases: ", x$cfrnow$n_cases,
          "   deaths by cut-off: ", x$cfrnow$n_deaths,
          "   naive CFR: ", round(attr(s, "naive_cfr"), 3))
  print(s)
  if (isTRUE(attr(s, "cfr_low_information"))) {
    message("  ! CFR only weakly identified (posterior close to prior); ",
            "few resolved deaths - interpret with caution.")
  }
  if (any(s$rhat > 1.01, na.rm = TRUE)) {
    message("  ! some Rhat > 1.01 - chains may not have converged.")
  }
  invisible(x)
}

#' Naive deaths / cases ratio
#'
#' @param n_deaths,n_cases Death and case counts.
#' @return The naive ratio, or `NA` when there are no cases.
#' @noRd
naive_cfr <- function(n_deaths, n_cases) {
  if (n_cases == 0) NA_real_ else n_deaths / n_cases
}

#' Mean and sd (days) of a delay from its brms native parameters
#'
#' `loc_var` is the location intercept and `scale_var` the second parameter's
#' intercept; the second parameter is log-linked. Handles lognormal (loc is
#' meanlog) and gamma (loc is log-mean, second parameter is shape).
#' @noRd
.delay_moments <- function(dr, loc_var, scale_var, family) {
  loc <- dr[[loc_var]]
  sc <- exp(dr[[scale_var]])
  if (family == "lognormal") {
    dmean <- exp(loc + sc^2 / 2)
    list(mean = dmean, sd = sqrt(exp(sc^2) - 1) * dmean)
  } else {
    dmean <- exp(loc)
    list(mean = dmean, sd = dmean / sqrt(sc))
  }
}

.cfr_quantities <- function(object) {
  dr <- posterior::as_draws_df(object)
  if (!"b_cfr_Intercept" %in% posterior::variables(dr)) {
    stop("summary() supports intercept-only `cfr` fits; for a `cfr ~ x` fit ",
         "use brms::posterior_epred() on the fit directly.", call. = FALSE)
  }
  fam <- object$cfrnow$family
  scale2 <- if (fam == "lognormal") "sigma" else "shape"
  d <- .delay_moments(dr, "b_Intercept",
                      paste0("b_", scale2, "_Intercept"), fam)
  res <- data.frame(cfr = stats::plogis(dr[["b_cfr_Intercept"]]),
                    delay_mean = d$mean, delay_sd = d$sd)
  if (isTRUE(object$cfrnow$use_recovery)) {
    rfam <- object$cfrnow$recovery_family
    rscale2 <- if (rfam == "lognormal") "sigma" else "shape"
    r <- .delay_moments(dr, "b_rmu_Intercept",
                        paste0("b_r", rscale2, "_Intercept"), rfam)
    res$recovery_mean <- r$mean
    res$recovery_sd <- r$sd
  }
  res$.chain <- dr$.chain
  res$.iteration <- dr$.iteration
  res$.draw <- dr$.draw
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
  qs <- .cfr_quantities(object)
  qcols <- paste0("q", probs * 100)
  sm <- posterior::summarise_draws(
    qs, mean = mean,
    stats::setNames(lapply(probs, function(p) {
      function(x) stats::quantile(x, p, names = FALSE)
    }), qcols),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  )
  out <- as.data.frame(sm)
  names(out)[1] <- "quantity"

  cfr_post_sd <- stats::sd(posterior::extract_variable(qs, "cfr"))
  prior_sd <- object$cfrnow$cfr_prior_sd
  low_info <- !is.na(prior_sd) && cfr_post_sd > info_tol * prior_sd
  attr(out, "naive_cfr") <- naive_cfr(object$cfrnow$n_deaths,
                                      object$cfrnow$n_cases)
  attr(out, "n_cases") <- object$cfrnow$n_cases
  attr(out, "n_deaths") <- object$cfrnow$n_deaths
  attr(out, "cfr_prior_sd") <- prior_sd
  attr(out, "cfr_low_information") <- low_info
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

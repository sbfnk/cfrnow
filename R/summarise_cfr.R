#' Naive deaths / cases ratio
#'
#' Single source of truth for the naive-CFR formula.
#' @param data A `cfrnow_data` list.
#' @return The naive ratio, or `NA` when there are no cases.
#' @noRd
naive_cfr <- function(data) {
  if (data$n_cases == 0) NA_real_ else data$n_deaths / data$n_cases
}

#' Prior sd of a Beta(a, b) distribution
#'
#' Used for the CFR low-information flag.
#' @param a,b Beta shape parameters.
#' @return The prior standard deviation.
#' @noRd
beta_sd <- function(a, b) {
  sqrt(a * b / ((a + b)^2 * (a + b + 1)))
}

#' Summarise a mixture-cure CFR fit
#'
#' Pools the posterior draws of the corrected CFR and the onset-to-death delay
#' (mean and standard deviation, in days) and reports quantile summaries with
#' convergence diagnostics (`rhat`, `ess_bulk`). The naive `deaths / cases`
#' ratio is returned as an attribute for comparison; in real time it
#' underestimates the corrected CFR because not every fatal case has died by the
#' cut-off. The `delay_mean`/`delay_sd` summaries (days) are recovered from the
#' native parameters as generated quantities, so they are the same whatever
#' family the delay used; with a fixed delay they are constant.
#'
#' When few deaths have resolved (a young outbreak), the CFR is only weakly
#' identified and its posterior stays close to the prior. This is reported via
#' the `cfr_low_information` attribute: `TRUE` when the CFR posterior sd exceeds
#' `info_tol` times the prior sd.
#'
#' @param object A `cfrnow_fit` from [fit_cfr()].
#' @param probs Quantiles to report.
#' @param info_tol Low-information threshold: flag when the CFR posterior sd is
#'   more than this fraction of the prior sd. Defaults to 0.9.
#' @param ... Unused.
#' @return A data frame with one row per summarised quantity (`cfr`,
#'   `delay_mean`, `delay_sd`), carrying `naive_cfr`, `n_cases`, `n_deaths`,
#'   `cfr_prior_sd` and `cfr_low_information` attributes.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(delay = dist.spec::LogNormal(2.4, 0.5))
#' fit <- fit_cfr(prepare_cfr_data(ll),
#'                delay = dist.spec::LogNormal(2.4, 0.5),
#'                cfr_prior = dist.spec::Beta(1, 1))
#' summary(fit)
#' }
#' @export
summary.cfrnow_fit <- function(object, probs = c(0.025, 0.5, 0.975),
                               info_tol = 0.9, ...) {
  if (!inherits(object, "cfrnow_fit")) {
    stop("`object` must come from fit_cfr().", call. = FALSE)
  }
  vars <- c("cfr", "delay_mean", "delay_sd")
  if (isTRUE(object$use_recovery)) {
    vars <- c(vars, "recovery_mean", "recovery_sd")
  }
  draws <- object$fit$draws(variables = vars)
  qcols <- paste0("q", probs * 100)

  qrow <- function(name) {
    m <- posterior::extract_variable_matrix(draws, name)  # iterations x chains
    x <- as.numeric(m)
    qv <- stats::quantile(x, probs = probs, names = FALSE)
    # A fixed delay makes delay_mean/delay_sd constant; rhat/ess are undefined
    # there, so report NA rather than a warning.
    constant <- stats::sd(x) == 0
    data.frame(quantity = name, mean = mean(x),
               stats::setNames(as.list(qv), qcols),
               rhat = if (constant) NA_real_ else posterior::rhat(m),
               ess_bulk = if (constant) NA_real_ else posterior::ess_bulk(m),
               check.names = FALSE)
  }
  out <- do.call(rbind, lapply(vars, qrow))

  cfr_draws <- as.numeric(posterior::extract_variable(draws, "cfr"))
  cfr_post_sd <- stats::sd(cfr_draws)
  cfr_prior_sd <- beta_sd(object$cfr_prior_shapes[["a"]],
                          object$cfr_prior_shapes[["b"]])

  attr(out, "naive_cfr") <- naive_cfr(object$data)
  attr(out, "n_cases") <- object$data$n_cases
  attr(out, "n_deaths") <- object$data$n_deaths
  attr(out, "cfr_prior_sd") <- cfr_prior_sd
  attr(out, "cfr_low_information") <- cfr_post_sd > info_tol * cfr_prior_sd
  out
}

#' @export
print.cfrnow_fit <- function(x, ...) {
  s <- summary(x)
  message("<cfrnow_fit> ", dist.spec::get_distribution(x$delay), " delay")
  message("  cases: ", x$data$n_cases,
          "   deaths by cut-off: ", x$data$n_deaths,
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

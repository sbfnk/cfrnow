#' Summarise a mixture-cure CFR fit
#'
#' Pools the posterior draws of the corrected CFR and the onset-to-death delay
#' (mean and standard deviation, in days) and reports quantile summaries. The
#' naive `deaths / cases` ratio is returned as an attribute for comparison; in
#' real time it underestimates the corrected CFR because not every fatal case
#' has died by the cut-off. The delay summaries are family-independent because
#' the model is parameterised directly by the delay mean and sd.
#'
#' @param object A `curecfr_fit` from [fit_cfr()].
#' @param probs Quantiles to report.
#' @param ... Unused.
#' @return A data frame with one row per summarised quantity (`cfr`,
#'   `delay_mean`, `delay_sd`) and `naive_cfr`, `n_cases`, `n_deaths`
#'   attributes.
#' @examples
#' \dontrun{
#' fit <- fit_cfr(prepare_cfr_data(simulate_linelist()))
#' summarise_cfr(fit)
#' }
#' @export
summarise_cfr <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  if (!inherits(object, "curecfr_fit")) {
    stop("`object` must come from fit_cfr().", call. = FALSE)
  }
  vars <- c("cfr", "delay_mean", "delay_sd")
  draws <- object$fit$draws(variables = vars, format = "draws_matrix")

  d <- object$data
  naive <- if (d$n_cases == 0) NA_real_ else d$n_deaths / d$n_cases

  qrow <- function(name) {
    x <- as.numeric(draws[, name])
    q <- stats::quantile(x, probs = probs, names = FALSE)
    data.frame(quantity = name, mean = mean(x),
               t(stats::setNames(q, paste0("q", probs * 100))),
               check.names = FALSE)
  }
  out <- do.call(rbind, lapply(vars, qrow))
  attr(out, "naive_cfr") <- naive
  attr(out, "n_cases") <- d$n_cases
  attr(out, "n_deaths") <- d$n_deaths
  out
}

#' @export
print.curecfr_fit <- function(x, ...) {
  naive <- if (x$data$n_cases == 0) NA_real_ else x$data$n_deaths / x$data$n_cases
  message("<curecfr_fit> ", x$delay_family, " delay")
  message("  cases: ", x$data$n_cases,
          "   deaths by cut-off: ", x$data$n_deaths,
          "   naive CFR: ", round(naive, 3))
  print(summarise_cfr(x))
  invisible(x)
}

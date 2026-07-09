#' Default CFR prior
#'
#' A weakly-informative logit-normal prior on the case fatality ratio, roughly
#' flat on \[0, 1]. Used by [fit_cfr()] when no `prior` is supplied.
#' @return A `brmsprior`.
#' @family fit
#' @export
cfr_default_prior <- function() {
  brms::set_prior("normal(0, 1.5)", class = "Intercept", dpar = "cfr")
}

#' Fit the real-time mixture-cure CFR model
#'
#' A thin wrapper over [epidist::epidist()] for the cfrnow cure model. Because
#' the model is an `epidist` subclass, both the CFR and the delay accept `brms`
#' formulas: pass a `brms::bf()` to put covariates on either, e.g.
#' `fit_cfr(d, bf(mu ~ 1, cfr ~ age))`.
#'
#' @param data Output of [prepare_cfr_data()], or an `epidist_cure_model` /
#'   data frame with `y`, `outcome`, `pwindow`, `swindow`.
#' @param formula A `stats::formula` or `brms::bf()`. A formula for `mu` (the
#'   delay location) is required; add a `cfr ~ ...` part to model the CFR.
#'   Defaults to `mu ~ 1`.
#' @param family The delay family: [brms::lognormal()] (default) or
#'   [brms::Gamma()].
#' @param prior `brms` priors. Defaults to [cfr_default_prior()]; supply your
#'   own (e.g. adding delay-parameter priors) to override.
#' @param ... Passed to [epidist::epidist()] and on to [brms::brm()]
#'   (e.g. `chains`, `iter`, `backend`, `seed`).
#' @return A `brmsfit` with class `cfrnow_fit`; summarise with [summary()].
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
#' d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
#' fit <- fit_cfr(d, backend = "cmdstanr")
#' summary(fit)
#' }
#' @family fit
#' @export
fit_cfr <- function(data, formula = mu ~ 1, family = brms::lognormal(),
                    prior = cfr_default_prior(), ...) {
  cure <- as_epidist_cure_model(data)
  fit <- epidist::epidist(cure, formula = formula, family = family,
                          prior = prior, merge_priors = FALSE, ...)
  fit$cfrnow <- list(
    n_cases = nrow(cure),
    n_deaths = sum(cure$outcome == .CURE_DEATH),
    family = family$family,
    use_recovery = isTRUE(attr(cure, "use_recovery")),
    cfr_prior_sd = .cfr_prior_sd(prior)
  )
  class(fit) <- c("cfrnow_fit", class(fit))
  fit
}

#' Prior sd of the CFR implied by a logit-scale Normal prior on its intercept
#'
#' Simulates the (0, 1)-scale sd from the `cfr` intercept prior so [summary()]
#' can flag when the posterior is barely tighter than the prior. Returns `NA`
#' when the prior is not a recognised `normal(mean, sd)` on the cfr intercept.
#' @param prior A `brmsprior`.
#' @return A numeric prior sd on the CFR scale, or `NA`.
#' @noRd
.cfr_prior_sd <- function(prior) {
  rows <- which(prior$dpar == "cfr" & prior$class == "Intercept")
  if (length(rows) != 1) return(NA_real_)
  m <- regmatches(
    prior$prior[rows],
    regexec("normal\\(([^,]+),\\s*([^)]+)\\)", prior$prior[rows])
  )[[1]]
  if (length(m) != 3) return(NA_real_)
  stats::sd(stats::plogis(stats::rnorm(1e5, as.numeric(m[2]), as.numeric(m[3]))))
}

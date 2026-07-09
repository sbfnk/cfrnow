#' CFR prior, specified in CFR units
#'
#' The CFR enters the model on the logit scale (it has covariates through the
#' formula), so its prior is a Normal on `logit(cfr)` rather than a `Beta`. This
#' helper lets you think in CFR units instead: it moment-matches a logit-Normal
#' whose induced mean and sd on the \[0, 1] CFR scale are approximately `mean`
#' and `sd`, and returns the corresponding `brms` prior on the `cfr` intercept.
#'
#' Rough guides: `cfr_prior(0.5, 0.29)` is near-flat (the default);
#' `cfr_prior(0.1, 0.1)` favours a low CFR; `cfr_prior(0.33, 0.12)` suits a
#' high-fatality pathogen. The CFR is weakly identified early in an outbreak, so
#' this prior can dominate — choose it deliberately.
#'
#' @param mean Prior mean of the CFR, on the \[0, 1] scale.
#' @param sd Prior sd of the CFR, on the \[0, 1] scale.
#' @return A `brmsprior` on the `cfr` intercept.
#' @family fit
#' @export
cfr_prior <- function(mean = 0.5, sd = 0.2) {
  moments <- function(m, s) {                          # of plogis(Normal(m, s))
    x <- seq(m - 6 * s, m + 6 * s, length.out = 2001)
    w <- stats::dnorm(x, m, s); w <- w / sum(w)
    p <- stats::plogis(x); mu <- sum(w * p)
    c(mu, sqrt(sum(w * (p - mu)^2)))
  }
  obj <- function(par) sum((moments(par[1], exp(par[2])) - c(mean, sd))^2)
  init <- c(stats::qlogis(mean), log(sd / (mean * (1 - mean) + 1e-6) + 1e-3))
  opt <- stats::optim(init, obj, method = "Nelder-Mead")
  brms::set_prior(sprintf("normal(%.4f, %.4f)", opt$par[1], exp(opt$par[2])),
                  class = "Intercept", dpar = "cfr")
}

#' Default CFR prior
#'
#' A weakly-informative CFR prior, roughly flat on \[0, 1]
#' (`cfr_prior(0.5, 0.29)`). Used by [fit_cfr()] when no `prior` is supplied.
#' @return A `brmsprior`.
#' @family fit
#' @export
cfr_default_prior <- function() {
  brms::set_prior("normal(0, 1.5)", class = "Intercept", dpar = "cfr")
}

#' Hold the onset-to-death delay fixed (the Ghani/Nishiura estimator)
#'
#' Returns `brms` `constant()` priors that pin the delay parameters, so only
#' `cfr` is estimated — the fixed-delay case fatality estimator. Supply the
#' delay by mean and sd (days); the family sets the parameterisation.
#'
#' @param family The delay family: [brms::lognormal()] or [brms::Gamma()].
#' @param mean,sd Mean and sd of the delay, in days.
#' @return A `brmsprior` (two `constant()` rows) for the delay parameters.
#' @family fit
#' @export
fix_delay <- function(family = brms::lognormal(), mean, sd) {
  fam <- brms:::validate_family(family)$family # nolint
  if (fam == "lognormal") {
    v <- log1p((sd / mean)^2)                          # sdlog^2
    loc <- log(mean) - v / 2                            # meanlog (identity link)
    scale <- log(sqrt(v))                               # log(sdlog)
    scale_dpar <- "sigma"
  } else if (fam == "gamma") {
    loc <- log(mean)                                    # log(mean)
    scale <- log((mean / sd)^2)                         # log(shape)
    scale_dpar <- "shape"
  } else {
    stop("cfrnow supports lognormal() and Gamma() delays only.", call. = FALSE)
  }
  c(brms::set_prior(sprintf("constant(%.8f)", loc), class = "Intercept"),
    brms::set_prior(sprintf("constant(%.8f)", scale), class = "Intercept",
                    dpar = scale_dpar))
}

#' Fit the real-time mixture-cure CFR model
#'
#' A thin wrapper over [epidist::epidist()] for the cfrnow cure model. Because
#' the model is an `epidist` subclass, both the CFR and the delay accept `brms`
#' formulas: pass a `brms::bf()` to put covariates on either, e.g.
#' `fit_cfr(d, bf(mu ~ 1, cfr ~ age))`.
#'
#' Define the CFR prior with [cfr_prior()] (in CFR units) or any `brms` prior on
#' the logit-scale `cfr` intercept. Add [fix_delay()] to `prior` to hold the
#' delay fixed. When the data carry recoveries, a two-outcome fit also estimates
#' an onset-to-recovery delay, whose family is `recovery_family` (defaulting to
#' `family`).
#'
#' @param data Output of [prepare_cfr_data()], or an `epidist_cure_model` /
#'   data frame with `y`, `outcome`, `pwindow`, `swindow`.
#' @param formula A `stats::formula` or `brms::bf()`. A formula for `mu` (the
#'   delay location) is required; add a `cfr ~ ...` part to model the CFR.
#'   Defaults to `mu ~ 1`.
#' @param family The delay family: [brms::lognormal()] (default) or
#'   [brms::Gamma()].
#' @param recovery_family Family for the onset-to-recovery delay in a
#'   two-outcome fit; defaults to `family`. Ignored without recovery data.
#' @param prior `brms` priors. Defaults to [cfr_default_prior()].
#' @param ... Passed to [epidist::epidist()] and on to [brms::brm()]
#'   (e.g. `chains`, `iter`, `backend`, `seed`).
#' @return A `brmsfit` with class `cfrnow_fit`; summarise with [summary()].
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
#' d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
#' fit <- fit_cfr(d, prior = cfr_prior(0.4, 0.15), backend = "cmdstanr")
#' summary(fit)
#' }
#' @family fit
#' @export
fit_cfr <- function(data, formula = mu ~ 1, family = brms::lognormal(),
                    recovery_family = family, prior = cfr_default_prior(), ...) {
  cure <- as_epidist_cure_model(data)
  if (isTRUE(attr(cure, "use_recovery"))) {
    attr(cure, "recovery_family") <- brms:::validate_family(recovery_family) # nolint
  }
  fit <- epidist::epidist(cure, formula = formula, family = family,
                          prior = prior, merge_priors = FALSE, ...)
  fit$cfrnow <- list(
    n_cases = nrow(cure),
    n_deaths = sum(cure$outcome == .CURE_DEATH),
    family = family$family,
    use_recovery = isTRUE(attr(cure, "use_recovery")),
    recovery_family = if (isTRUE(attr(cure, "use_recovery"))) {
      brms:::validate_family(recovery_family)$family # nolint
    } else NA_character_,
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

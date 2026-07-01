#' Priors for the mixture-cure CFR model
#'
#' Priors are placed on the interpretable onset-to-death delay **mean** and
#' **standard deviation** (in days) and on the CFR, so they do not depend on the
#' chosen delay family. Defaults are the BDBV/Ebolavirus values used by the
#' Julia reference model: an onset-to-death delay centred on the Isiro 2012
#' line-list reanalysis (mean ~= 12.75 d, sd ~= 7 d) and an EVD/BDBV
#' `Beta(6.6, 13.4)` CFR prior (mean ~= 0.33). The delay mean and sd carry
#' half-normal priors (truncated at zero).
#'
#' @param delay_mean_mean,delay_mean_sd Half-normal prior on the delay mean (days).
#' @param delay_sd_mean,delay_sd_sd Half-normal prior on the delay sd (days).
#' @param cfr_a,cfr_b Beta prior on the case fatality ratio.
#' @return A named list of prior hyperparameters for [fit_cfr()].
#' @examples
#' curecfr_priors(delay_mean_mean = 10)
#' @export
curecfr_priors <- function(delay_mean_mean = 12.75, delay_mean_sd = 3,
                           delay_sd_mean = 7, delay_sd_sd = 2,
                           cfr_a = 6.6, cfr_b = 13.4) {
  list(delay_mean_mean = delay_mean_mean, delay_mean_sd = delay_mean_sd,
       delay_sd_mean = delay_sd_mean, delay_sd_sd = delay_sd_sd,
       cfr_a = cfr_a, cfr_b = cfr_b)
}

#' Compile the mixture-cure CFR Stan model
#'
#' Compiles `inst/stan/curecfr.stan`, resolving the vendored primarycensored
#' functions via the package's Stan include path. The compiled model is cached
#' by cmdstanr, so repeated calls are cheap.
#'
#' @param ... Passed to [cmdstanr::cmdstan_model()].
#' @return A `CmdStanModel` object.
#' @examples
#' \dontrun{
#' model <- curecfr_model()
#' }
#' @export
curecfr_model <- function(...) {
  stan_file <- system.file("stan", "curecfr.stan", package = "curecfr")
  if (!nzchar(stan_file)) {
    stop("could not locate the packaged Stan model.", call. = FALSE)
  }
  cmdstanr::cmdstan_model(stan_file, include_paths = dirname(stan_file), ...)
}

#' Fit the real-time mixture-cure CFR model
#'
#' @param data A `curecfr_data` list from [prepare_cfr_data()].
#' @param delay_family Onset-to-death delay family; see [curecfr_families()].
#'   Defaults to `"gamma"` (the family selected for BDBV in the Isiro
#'   reanalysis). The mean/sd parameterisation makes priors and summaries
#'   identical across families.
#' @param priors Prior hyperparameters, see [curecfr_priors()].
#' @param model A compiled `CmdStanModel`; defaults to [curecfr_model()].
#' @param chains,parallel_chains,iter_warmup,iter_sampling,seed Passed to
#'   `CmdStanModel$sample()`.
#' @param init Initial values passed to `CmdStanModel$sample()`. Defaults to a
#'   generator drawing near the priors, which keeps warmup off the degenerate
#'   delay-sd = 0 boundary.
#' @param ... Further arguments to `CmdStanModel$sample()`.
#' @return A `curecfr_fit` object wrapping the `CmdStanMCMC` fit and `data`.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 300, cfr = 0.55)
#' d <- prepare_cfr_data(ll, obs_time = as.Date("2026-02-15"))
#' fit <- fit_cfr(d, delay_family = "gamma")
#' summarise_cfr(fit)
#' }
#' @export
fit_cfr <- function(data, delay_family = "gamma", priors = curecfr_priors(),
                    model = curecfr_model(),
                    chains = 4, parallel_chains = chains,
                    iter_warmup = 1000, iter_sampling = 1000,
                    seed = 20260508, init = NULL, ...) {
  if (!inherits(data, "curecfr_data")) {
    stop("`data` must come from prepare_cfr_data().", call. = FALSE)
  }
  dist_id <- delay_dist_id(delay_family)

  if (is.null(init)) {
    init <- function() {
      list(delay_mean = max(stats::rnorm(1, priors$delay_mean_mean,
                                         priors$delay_mean_sd / 2), 1),
           delay_sd = max(priors$delay_sd_mean, 1),
           cfr = stats::rbeta(1, priors$cfr_a, priors$cfr_b))
    }
  }

  stan_data <- c(
    data[c("n_death", "death_delay", "death_width",
           "n_cens", "censor_time", "censor_width", "n_resolved")],
    list(
      dist_id = dist_id,
      primary_id = primarycensored::pcd_stan_dist_id("uniform", type = "primary")
    ),
    priors
  )
  fit <- model$sample(
    data = stan_data, chains = chains, parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed, init = init, ...
  )
  structure(
    list(fit = fit, data = data, delay_family = delay_family),
    class = "curecfr_fit"
  )
}

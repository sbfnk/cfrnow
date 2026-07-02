#' Default onset-to-death delay specification
#'
#' A [dist.spec::LogNormal()] with `Normal()` priors on the native parameters,
#' reproducing the BDBV/Ebolavirus onset-to-death prior of the Julia reference
#' model (the Isiro 2012 line-list reanalysis): `meanlog ~ Normal(2.41, 0.2)`,
#' `sdlog ~ Normal(0.51, 0.15)` (delay mean ~= 12.75 d, sd ~= 7 d).
#'
#' Supply your own [dist.spec::LogNormal()] or [dist.spec::Gamma()] to
#' [fit_cfr()] to change the family or priors. Give each native parameter as a
#' `Normal()` to co-estimate it, or as a fixed number / [dist.spec::Fixed()] to
#' hold it fixed (fixing the whole delay gives the Ghani/Nishiura estimator).
#'
#' @return A `dist_spec` for the onset-to-death delay.
#' @examples
#' default_delay()
#' @export
default_delay <- function() {
  dist.spec::LogNormal(
    meanlog = dist.spec::Normal(mean = 2.41, sd = 0.2),
    sdlog = dist.spec::Normal(mean = 0.51, sd = 0.15)
  )
}

# Native parameter order per family (p1, p2), matching the Stan model.
delay_native_order <- function(family) {
  switch(
    family,
    lognormal = c("meanlog", "sdlog"),
    gamma = c("shape", "rate"),
    stop("unsupported delay family '", family,
         "'; cfrnow supports lognormal and gamma.", call. = FALSE)
  )
}

# Parse one native parameter of a delay `dist_spec` into Stan inputs: a fixed
# value (bare number or Fixed()) or a Normal() prior. dist.spec allows only
# Normal priors on parameters, which is what the Stan model expects.
parse_delay_param <- function(p, name) {
  if (is.numeric(p)) {
    return(list(est = 0L, fixed = p, prior_mean = 1, prior_sd = 1))
  }
  if (inherits(p, "dist_spec")) {
    dname <- dist.spec::get_distribution(p)
    pars <- dist.spec::get_parameters(p)
    if (dname == "fixed") {
      return(list(est = 0L, fixed = pars$value, prior_mean = 1, prior_sd = 1))
    }
    if (dname == "normal") {
      return(list(est = 1L, fixed = 1,
                  prior_mean = pars$mean, prior_sd = pars$sd))
    }
    stop("delay parameter '", name, "' must be fixed or a Normal() prior; got ",
         dname, ".", call. = FALSE)
  }
  stop("unrecognised specification for delay parameter '", name, "'.",
       call. = FALSE)
}

# Turn a delay `dist_spec` into the Stan data fields (dist_id + p1/p2 flags,
# fixed values and prior hyperparameters).
delay_to_stan_data <- function(delay) {
  if (!inherits(delay, "dist_spec")) {
    stop("`delay` must be a dist.spec distribution, ",
         "e.g. LogNormal() or Gamma().", call. = FALSE)
  }
  fam <- dist.spec::get_distribution(delay)
  native <- delay_native_order(fam)
  pars <- dist.spec::get_parameters(delay)
  if (!all(native %in% names(pars))) {
    stop("`delay` must be given in native parameters (", native[1], ", ",
         native[2], ") for a ", fam, " distribution.", call. = FALSE)
  }
  p1 <- parse_delay_param(pars[[native[1]]], native[1])
  p2 <- parse_delay_param(pars[[native[2]]], native[2])
  list(
    dist_id = primarycensored::pcd_stan_dist_id(fam, type = "delay"),
    primary_id = primarycensored::pcd_stan_dist_id("uniform", type = "primary"),
    p1_est = p1$est, p2_est = p2$est,
    p1_fixed = p1$fixed, p2_fixed = p2$fixed,
    p1_prior_mean = p1$prior_mean, p1_prior_sd = p1$prior_sd,
    p2_prior_mean = p2$prior_mean, p2_prior_sd = p2$prior_sd
  )
}

#' Compile the mixture-cure CFR Stan model
#'
#' Compiles `inst/stan/cfrnow.stan`, resolving the vendored primarycensored
#' functions via the package's Stan include path. The compiled model is cached
#' by cmdstanr, so repeated calls are cheap.
#'
#' @param ... Passed to [cmdstanr::cmdstan_model()].
#' @return A `CmdStanModel` object.
#' @examples
#' \dontrun{
#' model <- cfrnow_model()
#' }
#' @export
cfrnow_model <- function(...) {
  stan_file <- system.file("stan", "cfrnow.stan", package = "cfrnow")
  if (!nzchar(stan_file)) {
    stop("could not locate the packaged Stan model.", call. = FALSE)
  }
  cmdstanr::cmdstan_model(stan_file, include_paths = dirname(stan_file), ...)
}

#' Fit the real-time mixture-cure CFR model
#'
#' @param data A `cfrnow_data` list from [prepare_cfr_data()].
#' @param delay Onset-to-death delay as a dist.spec distribution
#'   ([dist.spec::LogNormal()] or [dist.spec::Gamma()]). Each native parameter
#'   is co-estimated when given as a `Normal()` prior, or held fixed when given
#'   as a number / [dist.spec::Fixed()]; fixing the whole delay yields the
#'   Ghani/Nishiura fixed-delay estimator. Defaults to [default_delay()].
#' @param cfr_a,cfr_b Beta prior on the case fatality ratio (default
#'   `Beta(6.6, 13.4)`, mean ~= 0.33).
#' @param model A compiled `CmdStanModel`; defaults to [cfrnow_model()].
#' @param chains,parallel_chains,iter_warmup,iter_sampling,seed,... Passed to
#'   `CmdStanModel$sample()`.
#' @return A `cfrnow_fit` object wrapping the `CmdStanMCMC` fit, `data` and the
#'   `delay` specification.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 300, cfr = 0.55)
#' d <- prepare_cfr_data(ll, obs_time = as.Date("2026-02-15"))
#' fit <- fit_cfr(d, delay = default_delay())
#' summarise_cfr(fit)
#' }
#' @export
fit_cfr <- function(data, delay = default_delay(), cfr_a = 6.6, cfr_b = 13.4,
                    model = cfrnow_model(),
                    chains = 4, parallel_chains = chains,
                    iter_warmup = 1000, iter_sampling = 1000,
                    seed = 20260508, ...) {
  if (!inherits(data, "cfrnow_data")) {
    stop("`data` must come from prepare_cfr_data().", call. = FALSE)
  }
  ds <- delay_to_stan_data(delay)

  stan_data <- c(
    data[c("n_death", "death_delay", "death_width",
           "n_cens", "censor_time", "censor_width", "n_resolved")],
    ds,
    list(cfr_a = cfr_a, cfr_b = cfr_b)
  )
  fit <- model$sample(
    data = stan_data, chains = chains, parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed, ...
  )
  structure(
    list(fit = fit, data = data, delay = delay,
         cfr_prior = c(a = cfr_a, b = cfr_b)),
    class = "cfrnow_fit"
  )
}

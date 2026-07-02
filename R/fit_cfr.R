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

# Turn a delay `dist_spec` into named Stan data fields: the family id (under
# `dist_id_name`) and the two native parameters under prefix `pfx` (fixed values
# and Normal-prior hyperparameters). Used for both the death (`p`) and recovery
# (`q`) delays.
stan_delay_fields <- function(delay, dist_id_name, pfx) {
  if (!inherits(delay, "dist_spec")) {
    stop("delay must be a dist.spec distribution, ",
         "e.g. LogNormal() or Gamma().", call. = FALSE)
  }
  fam <- dist.spec::get_distribution(delay)
  native <- delay_native_order(fam)
  pars <- dist.spec::get_parameters(delay)
  if (!all(native %in% names(pars))) {
    stop("delay must be given in native parameters (", native[1], ", ",
         native[2], ") for a ", fam, " distribution.", call. = FALSE)
  }
  p1 <- parse_delay_param(pars[[native[1]]], native[1])
  p2 <- parse_delay_param(pars[[native[2]]], native[2])
  out <- stats::setNames(
    list(primarycensored::pcd_stan_dist_id(fam, type = "delay"),
         p1$est, p2$est, p1$fixed, p2$fixed,
         p1$prior_mean, p1$prior_sd, p2$prior_mean, p2$prior_sd),
    c(dist_id_name,
      paste0(pfx, c("1_est", "2_est", "1_fixed", "2_fixed", "1_prior_mean",
                    "1_prior_sd", "2_prior_mean", "2_prior_sd")))
  )
  out
}

# A placeholder delay for the recovery slot when recovery is not modelled; its
# fixed values are never read (the Stan model gates them on `use_recovery`).
dummy_delay <- function() dist.spec::LogNormal(meanlog = 1, sdlog = 1)

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
#'   ([dist.spec::LogNormal()] or [dist.spec::Gamma()]); required, with no
#'   default, because a sensible delay is pathogen- and setting-specific. Each
#'   native parameter is co-estimated when given as a `Normal()` prior, or held
#'   fixed when given as a number / [dist.spec::Fixed()]; fixing the whole delay
#'   yields the Ghani/Nishiura fixed-delay estimator. For example the BDBV/Isiro
#'   onset-to-death prior is a `LogNormal` with `meanlog ~ Normal(2.41, 0.2)`,
#'   `sdlog ~ Normal(0.51, 0.15)` (mean ~= 12.75 d).
#' @param recovery_delay Optional onset-to-recovery delay (a dist.spec
#'   distribution, same conventions as `delay`). When supplied, cfrnow fits the
#'   competing-risks model that uses recovery *timing*: a recovered case
#'   contributes `(1 - cfr) * f_R(r)` and an unresolved case contributes
#'   `cfr * (1 - F_D(t)) + (1 - cfr) * (1 - F_R(t))`. When `NULL` (default),
#'   recovery timing is ignored and recovered cases contribute only the cure
#'   factor `1 - cfr`.
#' @param cfr_a,cfr_b Beta prior on the case fatality ratio (default
#'   `Beta(6.6, 13.4)`, mean ~= 0.33).
#' @param model A compiled `CmdStanModel`; defaults to [cfrnow_model()].
#' @param chains,parallel_chains,iter_warmup,iter_sampling,seed,... Passed to
#'   `CmdStanModel$sample()`.
#' @return A `cfrnow_fit` object wrapping the `CmdStanMCMC` fit, `data`, the
#'   `delay`/`recovery_delay` specifications and whether recovery was modelled.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 300, cfr = 0.55, recovery = TRUE)
#' d <- prepare_cfr_data(ll, obs_time = as.Date("2026-02-15"))
#' onset_to_death <- dist.spec::LogNormal(
#'   meanlog = dist.spec::Normal(2.41, 0.2),
#'   sdlog = dist.spec::Normal(0.51, 0.15)
#' )
#' fit <- fit_cfr(d, delay = onset_to_death)
#' summarise_cfr(fit)
#' }
#' @export
fit_cfr <- function(data, delay, recovery_delay = NULL,
                    cfr_a = 6.6, cfr_b = 13.4, model = cfrnow_model(),
                    chains = 4, parallel_chains = chains,
                    iter_warmup = 1000, iter_sampling = 1000,
                    seed = 20260508, ...) {
  if (!inherits(data, "cfrnow_data")) {
    stop("`data` must come from prepare_cfr_data().", call. = FALSE)
  }
  if (missing(delay)) {
    stop("supply a `delay` (an onset-to-death dist.spec distribution); ",
         "there is no default because it is pathogen- and setting-specific.",
         call. = FALSE)
  }
  use_recovery <- !is.null(recovery_delay)

  stan_data <- c(
    data[c("n_death", "death_delay", "death_width",
           "n_recovery", "recovery_delay", "recovery_width",
           "n_cens", "censor_time", "censor_width", "n_resolved")],
    stan_delay_fields(delay, "dist_id", "p"),
    stan_delay_fields(if (use_recovery) recovery_delay else dummy_delay(),
                      "recovery_dist_id", "q"),
    list(
      primary_id = primarycensored::pcd_stan_dist_id("uniform",
                                                     type = "primary"),
      use_recovery = as.integer(use_recovery),
      cfr_a = cfr_a, cfr_b = cfr_b
    )
  )
  fit <- model$sample(
    data = stan_data, chains = chains, parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed, ...
  )
  structure(
    list(fit = fit, data = data, delay = delay,
         recovery_delay = recovery_delay, use_recovery = use_recovery,
         cfr_prior = c(a = cfr_a, b = cfr_b)),
    class = "cfrnow_fit"
  )
}

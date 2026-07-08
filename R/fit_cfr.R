#' Native parameter order for a delay family
#'
#' The two native parameters in the positional order primarycensored's Stan
#' functions expect (p1, p2). Doubles as cfrnow's supported-family guard.
#' @param family Delay family name.
#' @return Character vector of the two native parameter names.
#' @noRd
delay_native_order <- function(family) {
  switch(
    family,
    lognormal = c("meanlog", "sdlog"),
    gamma = c("shape", "rate"),
    stop("unsupported delay family '", family,
         "'; cfrnow supports lognormal and gamma.", call. = FALSE)
  )
}

#' Parse one native delay parameter into Stan inputs
#'
#' Turns a native parameter (a bare number, `Fixed()`, or a `Normal()` prior)
#' into the fixed value or prior hyperparameters the Stan model reads. distspec
#' allows only Normal priors on parameters, which is what the model expects.
#' @param p The parameter: numeric, or a `dist_spec` (`Fixed()` / `Normal()`).
#' @param name Parameter name, used in error messages.
#' @return A list with `est`, `fixed`, `prior_mean` and `prior_sd`.
#' @noRd
parse_delay_param <- function(p, name) {
  if (is.numeric(p)) {
    return(list(est = 0L, fixed = p, prior_mean = 1, prior_sd = 1))
  }
  if (inherits(p, "dist_spec")) {
    dname <- distspec::get_distribution(p)
    pars <- distspec::get_parameters(p)
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

#' Turn a delay `dist_spec` into named Stan data fields
#'
#' The family id (under `dist_id_name`) and the two native parameters under
#' prefix `pfx` (fixed values and Normal-prior hyperparameters). Used for both
#' the death (`p`) and recovery (`q`) delays.
#' @param delay A distspec delay distribution.
#' @param dist_id_name Name of the family-id Stan field.
#' @param pfx Prefix for the two parameter fields.
#' @return A named list of Stan data fields.
#' @noRd
stan_delay_fields <- function(delay, dist_id_name, pfx) {
  if (!inherits(delay, "dist_spec")) {
    stop("delay must be a distspec distribution, ",
         "e.g. LogNormal() or Gamma().", call. = FALSE)
  }
  fam <- distspec::get_distribution(delay)
  native <- delay_native_order(fam)
  pars <- distspec::get_parameters(delay)
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

#' Placeholder recovery delay when recovery is not modelled
#'
#' Its fixed values are never read; the Stan model gates them on `use_recovery`.
#' @return A distspec LogNormal.
#' @noRd
dummy_delay <- function() distspec::LogNormal(meanlog = 1, sdlog = 1)

#' Default initial values that keep the sampler off degenerate delays
#'
#' cmdstan's random inits can start a delay at a near-zero mean, where the
#' survival probability `1 - F(t)` of every censored case is ~0 and the joint
#' log density is `-Inf` (cmdstan then warns "Error evaluating the log
#' probability at the initial value"). This starts each *estimated* native
#' parameter at its prior mean and the CFR at its prior mean, jittered per
#' chain, all on the natural scale; fixed parameters carry no sampled value and
#' need none.
#' @param stan_data The assembled Stan data list.
#' @return A function suitable for cmdstanr's `init` argument.
#' @noRd
cfr_stan_init <- function(stan_data) {
  function(chain_id = 1) {
    jit <- function(m) as.array(max(m, 0.01) * stats::runif(1, 0.9, 1.1))
    cfr0 <- stan_data$cfr_a / (stan_data$cfr_a + stan_data$cfr_b)
    init <- list(
      cfr = min(max(cfr0 * stats::runif(1, 0.9, 1.1), 1e-3), 1 - 1e-3)
    )
    if (stan_data$p1_est == 1) init$p1_par <- jit(stan_data$p1_prior_mean)
    if (stan_data$p2_est == 1) init$p2_par <- jit(stan_data$p2_prior_mean)
    if (stan_data$use_recovery == 1 && stan_data$q1_est == 1) {
      init$q1_par <- jit(stan_data$q1_prior_mean)
    }
    if (stan_data$use_recovery == 1 && stan_data$q2_est == 1) {
      init$q2_par <- jit(stan_data$q2_prior_mean)
    }
    init
  }
}

#' Read the Beta shape parameters from a CFR prior
#'
#' The CFR prior is a `distspec::Beta()` with fixed (numeric) shape parameters;
#' returns them as the `c(a, b)` the Stan model reads. Rejects anything that is
#' not a fixed Beta so a mis-specified prior fails loudly rather than silently.
#' @param cfr_prior A distspec Beta distribution.
#' @return A named numeric vector `c(a, b)`.
#' @noRd
parse_cfr_prior <- function(cfr_prior) {
  if (!inherits(cfr_prior, "dist_spec")) {
    stop("`cfr_prior` must be a distspec distribution, ",
         "e.g. distspec::Beta(shape1 = 1, shape2 = 1).", call. = FALSE)
  }
  if (distspec::get_distribution(cfr_prior) != "beta") {
    stop("`cfr_prior` must be a Beta() distribution.", call. = FALSE)
  }
  pars <- distspec::get_parameters(cfr_prior)
  if (!all(vapply(pars[c("shape1", "shape2")], is.numeric, logical(1)))) {
    stop("`cfr_prior` must have fixed (numeric) shape parameters, not priors.",
         call. = FALSE)
  }
  c(a = pars$shape1, b = pars$shape2)
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
#' @param delay Onset-to-death delay as a distspec distribution
#'   ([distspec::LogNormal()] or [distspec::Gamma()]); required, with no
#'   default. Each native parameter is co-estimated when given as a `Normal()`
#'   prior, or held fixed when given as a number / [distspec::Fixed()]; fixing
#'   the whole delay yields the Ghani/Nishiura fixed-delay estimator. For
#'   example the BDBV/Isiro onset-to-death prior is a `LogNormal` with
#'   `meanlog ~ Normal(2.41, 0.2)` and `sdlog ~ Normal(0.51, 0.15)`
#'   (mean ~= 12.75 d).
#' @param recovery_delay Optional onset-to-recovery delay (a distspec
#'   distribution, same conventions as `delay`). When supplied, cfrnow fits the
#'   two-outcome mixture-cure model that also uses recovery *timing*: a fatal
#'   case (probability `cfr`) dies at `F_D` and a non-fatal case recovers at
#'   `F_R`, so a recovered case contributes `(1 - cfr) * f_R(r)` and an
#'   unresolved case contributes
#'   `cfr * (1 - F_D(t)) + (1 - cfr) * (1 - F_R(t))`.
#'   This is a mixture over the two outcomes, not a cause-specific-hazards
#'   competing-risks model. When `NULL` (default), recovery timing is ignored
#'   and recovered cases contribute only the cure factor `1 - cfr`. Note this
#'   mode assumes recoveries among the unresolved are recorded: an actually
#'   recovered case whose recovery is missing stays censored and, as time
#'   passes, is pushed toward the fatal branch, biasing `cfr` *up* (the mirror
#'   of incomplete death ascertainment). The death-only default is insensitive
#'   to this.
#' @param cfr_prior Prior on the case fatality ratio, as a
#'   [distspec::Beta()] with fixed shape parameters; required, with no default.
#'   Because the CFR is only weakly identified early in an outbreak (few deaths
#'   resolved), this prior can dominate, so choose it deliberately rather than
#'   reaching for a default. Some reference choices and what they imply (the
#'   Beta "prior sample size" is `shape1 + shape2`, so larger sums pull harder):
#'   `Beta(1, 1)` is uniform on \[0, 1] (weakly informative, mean 0.5);
#'   `Beta(1, 9)` has mean 0.1 and favours a low CFR;
#'   `Beta(6.6, 13.4)` has mean 0.33 and is appropriate only for a
#'   high-fatality pathogen (it is roughly the BDBV/Isiro prior). Specify by
#'   mean and sd instead if that is easier, e.g. `Beta(mean = 0.1, sd = 0.1)`.
#' @param model A compiled `CmdStanModel`; defaults to [cfrnow_model()].
#' @param chains,parallel_chains,iter_warmup,iter_sampling,seed,... Passed to
#'   `CmdStanModel$sample()`. `seed` defaults to `NULL`, so cmdstanr draws a
#'   random seed; set it for reproducible fits. Unless you pass your own `init`,
#'   the sampler starts each estimated delay parameter and the CFR at their
#'   prior means (jittered per chain), which avoids the degenerate
#'   near-zero-delay starts that make cmdstan reject the initial log density.
#' @return A `cfrnow_fit` object wrapping the `CmdStanMCMC` fit, `data`, the
#'   `delay`/`recovery_delay` specifications and whether recovery was modelled.
#' @examples
#' \dontrun{
#' ll <- simulate_linelist(n = 300, cfr = 0.55,
#'                         delay = distspec::LogNormal(mean = 12.75, sd = 7),
#'                         recovery = distspec::LogNormal(mean = 21, sd = 9))
#' d <- prepare_cfr_data(ll, obs_time = as.Date("2026-02-15"))
#' onset_to_death <- distspec::LogNormal(
#'   meanlog = distspec::Normal(2.41, 0.2),
#'   sdlog = distspec::Normal(0.51, 0.15)
#' )
#' fit <- fit_cfr(d, delay = onset_to_death,
#'                cfr_prior = distspec::Beta(6.6, 13.4))
#' summary(fit)
#' }
#' @export
fit_cfr <- function(data, delay, recovery_delay = NULL, cfr_prior,
                    model = cfrnow_model(),
                    chains = 4, parallel_chains = chains,
                    iter_warmup = 1000, iter_sampling = 1000,
                    seed = NULL, ...) {
  if (!inherits(data, "cfrnow_data")) {
    stop("`data` must come from prepare_cfr_data().", call. = FALSE)
  }
  if (missing(delay)) {
    stop("supply a `delay` (an onset-to-death distspec distribution).",
         call. = FALSE)
  }
  if (missing(cfr_prior)) {
    stop("supply a `cfr_prior` (a distspec Beta() on the CFR); ",
         "there is no default. See ?fit_cfr for reference choices.",
         call. = FALSE)
  }
  cfr_shapes <- parse_cfr_prior(cfr_prior)
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
      cfr_a = cfr_shapes[["a"]], cfr_b = cfr_shapes[["b"]]
    )
  )
  dots <- list(...)
  if (!("init" %in% names(dots))) dots$init <- cfr_stan_init(stan_data)
  fit <- do.call(model$sample, c(
    list(data = stan_data, chains = chains, parallel_chains = parallel_chains,
         iter_warmup = iter_warmup, iter_sampling = iter_sampling, seed = seed),
    dots
  ))
  structure(
    list(fit = fit, data = data, delay = delay,
         recovery_delay = recovery_delay, use_recovery = use_recovery,
         cfr_prior = cfr_prior,
         cfr_prior_shapes = c(a = cfr_shapes[["a"]], b = cfr_shapes[["b"]])),
    class = "cfrnow_fit"
  )
}

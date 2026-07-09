# Translate distspec delay / cfr priors into the brms family and priors the
# epidist engine reads. distspec is cfrnow's user-facing idiom; this keeps the
# brms parameterisation an implementation detail.

`%||%` <- function(a, b) if (is.null(a)) b else a # nolint: coalesce_linter.

# A native delay parameter is either a fixed number or a Normal() prior.
.param_spec <- function(p) {
  if (is.numeric(p)) {
    return(list(fixed = p))
  }
  if (inherits(p, "dist_spec")) {
    dname <- get_distribution(p)
    pr <- get_parameters(p)
    if (dname == "fixed") {
      return(list(fixed = pr$value))
    }
    if (dname == "normal") {
      return(list(mean = pr$mean, sd = pr$sd))
    }
  }
  stop("delay parameters must be fixed numbers or Normal() priors.",
       call. = FALSE)
}

# One brms prior row for a dpar. `log_link` transforms to the link scale: a
# fixed value becomes a constant() prior, a Normal() becomes a Normal on the
# link scale (delta method for the sd when log-linked).
.prior_row <- function(spec, dpar, log_link) {
  tf <- function(x) if (log_link) log(x) else x
  prior_args <- if (nzchar(dpar)) {
    list(class = "Intercept", dpar = dpar)
  } else {
    list(class = "Intercept")
  }
  code <- if (!is.null(spec$fixed)) {
    sprintf("constant(%.8f)", tf(spec$fixed))
  } else {
    sc <- if (log_link) spec$sd / spec$mean else spec$sd
    sprintf("normal(%.6f, %.6f)", tf(spec$mean), sc)
  }
  do.call(brms::set_prior, c(list(code), prior_args))
}

# distspec Gamma() uses (shape, rate); brms Gamma uses (mu = mean, shape), both
# log-linked. Build the mu prior from mean = shape / rate.
.gamma_mu_prior <- function(sh, rate_spec, dpar) {
  if (!is.null(sh$fixed) && !is.null(rate_spec$fixed)) {
    return(.prior_row(list(fixed = sh$fixed / rate_spec$fixed), dpar, TRUE))
  }
  shc <- sh$fixed %||% sh$mean
  rtc <- rate_spec$fixed %||% rate_spec$mean
  vlog <- sqrt(((sh$sd %||% 0) / shc)^2 + ((rate_spec$sd %||% 0) / rtc)^2)
  prior_args <- if (nzchar(dpar)) {
    list(class = "Intercept", dpar = dpar)
  } else {
    list(class = "Intercept")
  }
  do.call(brms::set_prior,
          c(list(sprintf("normal(%.6f, %.6f)", log(shc / rtc), vlog)),
            prior_args))
}

# A distspec delay -> list(family, prior). `main = TRUE` targets the death delay
# (mu is the main dpar); `main = FALSE` the recovery delay (r-prefixed dpars).
.delay_family_prior <- function(delay, main = TRUE) {
  if (!inherits(delay, "dist_spec")) {
    stop("`delay` must be a distspec distribution, ",
         "e.g. LogNormal() or Gamma().", call. = FALSE)
  }
  fam <- get_distribution(delay)
  pars <- get_parameters(delay)
  loc_dpar <- if (main) "" else "rmu"
  if (fam == "lognormal") {
    scale_dpar <- if (main) "sigma" else "rsigma"
    prior <- c(.prior_row(.param_spec(pars$meanlog), loc_dpar, FALSE),
               .prior_row(.param_spec(pars$sdlog), scale_dpar, TRUE))
    list(family = brms::lognormal(), prior = prior)
  } else if (fam == "gamma") {
    scale_dpar <- if (main) "shape" else "rshape"
    sh <- .param_spec(pars$shape)
    rate_spec <- .param_spec(pars$rate)
    prior <- c(.gamma_mu_prior(sh, rate_spec, loc_dpar),
               .prior_row(sh, scale_dpar, TRUE))
    list(family = stats::Gamma(link = "log"), prior = prior)
  } else {
    stop("cfrnow supports LogNormal() and Gamma() delays only.", call. = FALSE)
  }
}

# A Normal(m, s) on logit(cfr) whose induced mean/sd on the [0, 1] CFR scale
# match `mean`/`sd` (moment-matched).
.cfr_logitnormal <- function(cfr_mean, cfr_sd) {
  moments <- function(m, s) {
    x <- seq(m - 6 * s, m + 6 * s, length.out = 2001)
    w <- stats::dnorm(x, m, s)
    w <- w / sum(w)
    p <- stats::plogis(x)
    pm <- sum(w * p)
    c(pm, sqrt(sum(w * (p - pm)^2)))
  }
  obj <- function(par) {
    sum((moments(par[1], exp(par[2])) - c(cfr_mean, cfr_sd))^2)
  }
  init <- c(stats::qlogis(cfr_mean),
            log(cfr_sd / (cfr_mean * (1 - cfr_mean) + 1e-6) + 1e-3))
  opt <- stats::optim(init, obj, method = "Nelder-Mead")
  brms::set_prior(sprintf("normal(%.4f, %.4f)", opt$par[1], exp(opt$par[2])),
                  class = "Intercept", dpar = "cfr")
}

# A distspec Beta() CFR prior -> the matching logit-scale brms prior.
.cfr_prior_to_brms <- function(cfr_prior) {
  ok <- inherits(cfr_prior, "dist_spec") &&
    get_distribution(cfr_prior) == "beta"
  if (!ok) {
    stop("`cfr_prior` must be a distspec Beta() distribution.", call. = FALSE)
  }
  p <- get_parameters(cfr_prior)
  a <- p$shape1
  b <- p$shape2
  .cfr_logitnormal(a / (a + b), sqrt(a * b / ((a + b)^2 * (a + b + 1))))
}

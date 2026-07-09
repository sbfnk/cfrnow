#' Real-time CFR as an `epidist` mixture-cure model
#'
#' `cfrnow` registers a mixture-cure survival model as an [epidist::epidist()]
#' model type. Each case is fatal with probability `cfr` and, when fatal, dies
#' at an onset-to-death delay; cases still alive at the observation cut-off are
#' right-censored, which corrects the downward bias of the naive deaths / cases
#' ratio in real time. Because the model is an `epidist` subclass, the CFR and
#' the delay both take `brms` formulas, e.g.
#' `epidist(data, bf(mu ~ 1, cfr ~ age), family = lognormal())`.
#'
#' The delay distribution's location is `mu` (as `epidist` expects); the cure
#' probability `cfr` is an additional dpar with a logit link. Supported delay
#' families are `lognormal()` and `Gamma()`.
#'
#' @name cfrnow-cure-model
NULL

# Outcome codes carried in the `outcome` vreal (match cfrnow's likelihood):
#   1 = observed death (timed)      -> log(cfr) + delay density
#   3 = resolved non-death (untimed)-> log1m(cfr)  (cure factor)
#   0 = censored survivor           -> log1m(cfr * F_delay(t))
.CURE_DEATH <- 1
.CURE_RESOLVED <- 3
.CURE_CENSORED <- 0

#' Build an `epidist_cure_model` from a line list or `prepare_cfr_data()` output
#'
#' Produces the one-row-per-case data frame the model needs: `y` (the delay for
#' deaths, the follow-up time for survivors), an `outcome` code, and the primary
#' (`pwindow`) and secondary (`swindow`) censoring widths.
#'
#' @param data Either the list returned by [prepare_cfr_data()] or a data frame
#'   already carrying `y`, `outcome`, `pwindow` and `swindow`.
#' @return A data frame of class `epidist_cure_model`.
#' @family cure_model
#' @export
as_epidist_cure_model <- function(data) {
  if (inherits(data, "cfrnow_data")) {
    if (isTRUE(data$n_recovery > 0)) {
      stop("recovery-timing (two-outcome) fits are not yet available on the ",
           "epidist backend; refit without `recovery_date` / recovery_delay.",
           call. = FALSE)
    }
    deaths <- if (data$n_deaths > 0) {
      data.frame(y = as.integer(data$death_delay), outcome = .CURE_DEATH,
                 pwindow = data$death_width, swindow = 1)
    } else NULL
    resolved <- if (data$n_resolved > 0) {
      data.frame(y = 0L, outcome = .CURE_RESOLVED,
                 pwindow = 1, swindow = 1)[rep(1, data$n_resolved), ]
    } else NULL
    cens <- if (data$n_cens > 0) {
      data.frame(y = as.integer(data$censor_time), outcome = .CURE_CENSORED,
                 pwindow = data$censor_width, swindow = 1)
    } else NULL
    data <- rbind(deaths, resolved, cens)
  }
  stopifnot(all(c("y", "outcome", "pwindow", "swindow") %in% names(data)))
  data <- as.data.frame(data)
  class(data) <- c("epidist_cure_model", "data.frame")
  data
}

#' @method assert_epidist epidist_cure_model
#' @export
assert_epidist.epidist_cure_model <- function(data, ...) {
  checkmate::assert_data_frame(data)
  checkmate::assert_names(names(data),
    must.include = c("y", "outcome", "pwindow", "swindow"))
  checkmate::assert_integerish(data$y)
  checkmate::assert_subset(data$outcome, c(.CURE_DEATH, .CURE_RESOLVED,
                                           .CURE_CENSORED))
  invisible(TRUE)
}

#' @method epidist_family_model epidist_cure_model
#' @export
epidist_family_model.epidist_cure_model <- function(data, family, ...) {
  if (!family$family %in% c("lognormal", "gamma")) {
    stop("cfrnow supports lognormal() and Gamma() delays only.", call. = FALSE)
  }
  # Per-dpar links: mu from the family, other (positive) delay params default to
  # log when the base family did not carry an explicit link, then logit for cfr.
  other_dpars <- setdiff(family$dpars, "mu")
  other_links <- vapply(other_dpars, function(dp) {
    l <- family[[paste0("link_", dp)]]
    if (is.null(l) || !nzchar(l)) "log" else l
  }, character(1))
  brms::custom_family(
    paste0("cfrnow_", family$family),
    dpars = c(family$dpars, "cfr"),                 # delay dpars, then cfr
    links = c(family$link, other_links, "logit"),
    type = "int",
    vars = c("vreal1[n]", "vreal2[n]", "vreal3[n]", "primary_params"),
    loop = TRUE,
    log_lik = function(i, prep, ...) rep(NA_real_, prep$ndraws),
    posterior_predict = function(i, prep, ...) rep(NA_real_, prep$ndraws)
  )
}

#' @method epidist_formula_model epidist_cure_model
#' @export
epidist_formula_model.epidist_cure_model <- function(data, formula, ...) {
  stats::update(formula, y | vreal(outcome, pwindow, swindow) ~ .)
}

#' @method epidist_transform_data_model epidist_cure_model
#' @export
epidist_transform_data_model.epidist_cure_model <- function(data, family,
                                                            formula, ...) {
  data
}

#' @method epidist_model_prior epidist_cure_model
#' @export
epidist_model_prior.epidist_cure_model <- function(data, formula, ...) NULL

#' @method epidist_stancode epidist_cure_model
#' @export
epidist_stancode.epidist_cure_model <- function(data, family, formula, ...) {
  family_name <- sub("^cfrnow_", "", family$name)
  dist_id <- primarycensored::pcd_stan_dist_id(family_name, type = "delay")
  delay_dpars <- utils::head(family$dpars, -1)          # drop trailing "cfr"

  tmpl <- '
  real cfrnow_family_lpmf(data int y, dpars_A, real cfr, data real outcome,
                          data real pwindow, data real swindow,
                          array[] real primary_params) {
    if (outcome == 1) {
      return log(cfr) + primarycensored_lpmf(
          y | dist_id, {dpars_B}, pwindow, y + swindow, 0.0,
          positive_infinity(), primary_id, primary_params);
    } else if (outcome == 3) {
      return log1m(cfr);
    } else {
      real fbar = primarycensored_cdf(
          y | dist_id, {dpars_B}, pwindow, 0.0, positive_infinity(),
          primary_id, primary_params);
      return log1m(cfr * fbar);
    }
  }'
  tmpl <- gsub("family", family_name, tmpl, fixed = TRUE)
  tmpl <- gsub("dpars_A", toString(paste0("real ", delay_dpars)), tmpl,
               fixed = TRUE)
  tmpl <- gsub("dpars_B", family$param, tmpl, fixed = TRUE)
  tmpl <- gsub("dist_id", dist_id, tmpl, fixed = TRUE)
  tmpl <- gsub("primary_id",
               primarycensored::pcd_stan_dist_id("uniform", type = "primary"),
               tmpl, fixed = TRUE)

  brms::stanvar(scode = primarycensored::pcd_load_stan_functions(),
                block = "functions") +
    brms::stanvar(scode = tmpl, block = "functions") +
    brms::stanvar(scode = "array[0] real primary_params;", block = "parameters")
}

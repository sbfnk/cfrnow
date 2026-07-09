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
#' When the line list carries recovery dates, `prepare_cfr_data()` times those
#' recoveries and the fit becomes a two-outcome mixture-cure model that also uses
#' onset-to-recovery timing: a non-fatal case recovers at a second delay, so a
#' recovered case contributes `(1 - cfr) f_R(r)` and an unresolved case
#' `cfr (1 - F_D(t)) + (1 - cfr)(1 - F_R(t))`. The recovery delay shares the
#' death delay's family.
#'
#' The delay distribution's location is `mu` (as `epidist` expects); the cure
#' probability `cfr` is an additional dpar with a logit link. Supported delay
#' families are `lognormal()` and `Gamma()`.
#'
#' @name cfrnow-cure-model
NULL

# Outcome codes carried in the `outcome` vreal (match cfrnow's likelihood):
#   1 = observed death (timed)       -> log(cfr) + death delay density
#   2 = recorded recovery (timed)    -> log1m(cfr) + recovery delay density
#   3 = resolved non-death (untimed) -> log1m(cfr)  (cure factor)
#   0 = censored survivor            -> death-only: log1m(cfr F_D(t));
#                                       two-outcome: log_sum_exp mixture
.CURE_DEATH <- 1
.CURE_RECOVERY <- 2
.CURE_RESOLVED <- 3
.CURE_CENSORED <- 0

#' Build an `epidist_cure_model` from a line list or `prepare_cfr_data()` output
#'
#' Produces the one-row-per-case data frame the model needs: `y` (a delay for
#' deaths and recoveries, the follow-up time for survivors), an `outcome` code,
#' and the primary (`pwindow`) and secondary (`swindow`) censoring widths. When
#' the input carries timed recoveries the result is flagged for the two-outcome
#' fit via a `use_recovery` attribute.
#'
#' @param data Either the list returned by [prepare_cfr_data()] or a data frame
#'   already carrying `y`, `outcome`, `pwindow` and `swindow`.
#' @return A data frame of class `epidist_cure_model`.
#' @family cure_model
#' @export
as_epidist_cure_model <- function(data) {
  if (inherits(data, "cfrnow_data")) {
    rows <- function(n, y, code, width) {
      if (n > 0) {
        data.frame(y = as.integer(y), outcome = code, pwindow = width,
                   swindow = 1)
      } else NULL
    }
    resolved <- if (data$n_resolved > 0) {
      data.frame(y = 0L, outcome = .CURE_RESOLVED, pwindow = 1,
                 swindow = 1)[rep(1, data$n_resolved), ]
    } else NULL
    data <- rbind(
      rows(data$n_deaths, data$death_delay, .CURE_DEATH, data$death_width),
      rows(data$n_recovery, data$recovery_delay, .CURE_RECOVERY,
           data$recovery_width),
      resolved,
      rows(data$n_cens, data$censor_time, .CURE_CENSORED, data$censor_width)
    )
  }
  stopifnot(all(c("y", "outcome", "pwindow", "swindow") %in% names(data)))
  data <- as.data.frame(data)
  class(data) <- c("epidist_cure_model", "data.frame")
  attr(data, "use_recovery") <- any(data$outcome == .CURE_RECOVERY)
  data
}

#' @method assert_epidist epidist_cure_model
#' @export
assert_epidist.epidist_cure_model <- function(data, ...) {
  checkmate::assert_data_frame(data)
  checkmate::assert_names(names(data),
    must.include = c("y", "outcome", "pwindow", "swindow"))
  checkmate::assert_integerish(data$y)
  checkmate::assert_subset(data$outcome,
    c(.CURE_DEATH, .CURE_RECOVERY, .CURE_RESOLVED, .CURE_CENSORED))
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
  delay_links <- c(family$link, other_links)
  dpars <- c(family$dpars, "cfr")
  links <- c(delay_links, "logit")
  if (isTRUE(attr(data, "use_recovery"))) {          # recovery shares the family
    dpars <- c(dpars, paste0("r", family$dpars))
    links <- c(links, delay_links)
  }
  brms::custom_family(
    paste0("cfrnow_", family$family),
    dpars = dpars, links = links, type = "int",
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
  primary_id <- primarycensored::pcd_stan_dist_id("uniform", type = "primary")
  cfr_pos <- match("cfr", family$dpars)
  delay_dpars <- family$dpars[seq_len(cfr_pos - 1)]
  use_recovery <- cfr_pos < length(family$dpars)

  # death-side param declarations and the family's Stan reparameterisation
  DA <- toString(paste0("real ", delay_dpars))
  DB <- family$param
  # recovery side: same reparameterisation with r-prefixed parameter names
  RB <- DB
  for (dp in delay_dpars) RB <- gsub(paste0("\\b", dp, "\\b"), paste0("r", dp), RB)
  RA <- toString(paste0("real r", delay_dpars))

  body <- if (use_recovery) '
    if (outcome == 1) {
      return log(cfr) + primarycensored_lpmf(
          y | @DID@, {@DB@}, pwindow, y + swindow, 0.0, positive_infinity(),
          @PID@, primary_params);
    } else if (outcome == 2) {
      return log1m(cfr) + primarycensored_lpmf(
          y | @DID@, {@RB@}, pwindow, y + swindow, 0.0, positive_infinity(),
          @PID@, primary_params);
    } else if (outcome == 3) {
      return log1m(cfr);
    } else {
      real fbar_d = primarycensored_cdf(y | @DID@, {@DB@}, pwindow, 0.0,
          positive_infinity(), @PID@, primary_params);
      real fbar_r = primarycensored_cdf(y | @DID@, {@RB@}, pwindow, 0.0,
          positive_infinity(), @PID@, primary_params);
      return log_sum_exp(log(cfr) + log1m(fbar_d), log1m(cfr) + log1m(fbar_r));
    }' else '
    if (outcome == 1) {
      return log(cfr) + primarycensored_lpmf(
          y | @DID@, {@DB@}, pwindow, y + swindow, 0.0, positive_infinity(),
          @PID@, primary_params);
    } else if (outcome == 3) {
      return log1m(cfr);
    } else {
      real fbar = primarycensored_cdf(y | @DID@, {@DB@}, pwindow, 0.0,
          positive_infinity(), @PID@, primary_params);
      return log1m(cfr * fbar);
    }'

  sig <- if (use_recovery) {
    sprintf("data int y, %s, real cfr, %s, data real outcome,\n%s", DA, RA,
            "  data real pwindow, data real swindow, array[] real primary_params")
  } else {
    sprintf("data int y, %s, real cfr, data real outcome,\n%s", DA,
            "  data real pwindow, data real swindow, array[] real primary_params")
  }
  tmpl <- sprintf("\n  real cfrnow_%s_lpmf(%s) {%s\n  }",
                  family_name, sig, body)
  for (r in list(c("@DID@", dist_id), c("@PID@", primary_id),
                 c("@DB@", DB), c("@RB@", RB))) {
    tmpl <- gsub(r[1], r[2], tmpl, fixed = TRUE)
  }

  brms::stanvar(scode = primarycensored::pcd_load_stan_functions(),
                block = "functions") +
    brms::stanvar(scode = tmpl, block = "functions") +
    brms::stanvar(scode = "array[0] real primary_params;", block = "parameters")
}

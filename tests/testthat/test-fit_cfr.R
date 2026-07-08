skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  ok <- tryCatch(!is.null(cmdstanr::cmdstan_version()), error = function(e) FALSE)
  testthat::skip_if(!isTRUE(ok), "cmdstan not installed")
}

test_that("retrospective fit recovers a known CFR", {
  skip_if_no_cmdstan()
  set.seed(1)
  ll <- simulate_linelist(n = 600, cfr = 0.45,
                          delay = Gamma(mean = 12.75, sd = 7))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  delay <- LogNormal(meanlog = Normal(2.41, 0.2),
                                sdlog = Normal(0.51, 0.15))
  fit <- fit_cfr(d, delay = delay, cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0)
  s <- summary(fit)
  cfr_row <- s[s$quantity == "cfr", ]
  expect_gt(cfr_row[["q97.5"]], 0.45)
  expect_lt(cfr_row[["q2.5"]], 0.45)
})

test_that("summary carries convergence diagnostics and identifiability flag", {
  skip_if_no_cmdstan()
  set.seed(3)
  ll <- simulate_linelist(n = 500, cfr = 0.5,
                          delay = Gamma(mean = 12.75, sd = 7))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  delay <- LogNormal(meanlog = Normal(2.41, 0.2),
                                sdlog = Normal(0.51, 0.15))
  fit <- fit_cfr(d, delay = delay, cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0)
  s <- summary(fit)
  expect_true(all(c("rhat", "ess_bulk") %in% names(s)))
  expect_true(all(s$rhat < 1.05))                       # well-mixed on ample data
  expect_false(is.null(attr(s, "cfr_low_information")))
  expect_false(attr(s, "cfr_low_information"))           # 500 cases: informative
})

test_that("a fixed delay runs the fixed-F (Ghani/Nishiura) estimator", {
  skip_if_no_cmdstan()
  set.seed(4)
  ll <- simulate_linelist(n = 400, cfr = 0.5,
                          delay = Gamma(mean = 12.75, sd = 7))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  fit <- fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51),
                 cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0)
  s <- summary(fit)
  # a fixed delay is constant, so its diagnostics are NA; cfr is still estimated
  expect_true(is.na(s[s$quantity == "delay_mean", "rhat"]))
  expect_gt(s[s$quantity == "cfr", "q97.5"], 0.5)
  expect_lt(s[s$quantity == "cfr", "q2.5"], 0.5)
})

test_that("real-time correction lifts the estimate above the naive ratio", {
  skip_if_no_cmdstan()
  set.seed(2)
  # Growing epidemic sampled mid-outbreak: many recent, not-yet-resolved cases.
  ll <- simulate_linelist(n = 800, cfr = 0.6, onset_days = 40,
                          delay = Gamma(mean = 12.75, sd = 7))
  cut <- max(ll$onset_date) - 2          # cut-off soon after the last onsets
  d <- prepare_cfr_data(ll, obs_time = cut)
  naive <- d$n_deaths / d$n_cases
  delay <- LogNormal(meanlog = Normal(2.41, 0.2),
                                sdlog = Normal(0.51, 0.15))
  fit <- fit_cfr(d, delay = delay, cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0)
  s <- summary(fit)
  cfr_med <- s[s$quantity == "cfr", "q50"]
  expect_gt(cfr_med, naive)             # corrected > downward-biased naive
})

test_that("competing-risks fit uses recovery timing and recovers CFR + F_R", {
  skip_if_no_cmdstan()
  set.seed(5)
  ll <- simulate_linelist(n = 700, cfr = 0.5,
                          delay = Gamma(mean = 12.75, sd = 7),
                          recovery = Gamma(mean = 21, sd = 9))
  d <- prepare_cfr_data(ll, obs_time = NULL)   # retrospective, all resolved
  expect_gt(d$n_recovery, 0)
  fit <- fit_cfr(
    d,
    delay = Gamma(shape = Normal(3.3, 1),
                             rate = Normal(0.26, 0.08)),
    recovery_delay = Gamma(shape = Normal(5, 2),
                                      rate = Normal(0.24, 0.08)),
    cfr_prior = Beta(1, 1),
    chains = 2, parallel_chains = 2,
    iter_warmup = 500, iter_sampling = 500, refresh = 0
  )
  s <- summary(fit)
  expect_true(all(c("recovery_mean", "recovery_sd") %in% s$quantity))
  cfr <- s[s$quantity == "cfr", ]
  expect_gt(cfr[["q97.5"]], 0.5)
  expect_lt(cfr[["q2.5"]], 0.5)
  rmean <- s[s$quantity == "recovery_mean", ]  # true onset-to-recovery mean = 21
  expect_gt(rmean[["q97.5"]], 21)
  expect_lt(rmean[["q2.5"]], 21)
})

test_that("print.cfrnow_fit reports the delay, counts and a summary", {
  skip_if_no_cmdstan()
  set.seed(6)
  ll <- simulate_linelist(n = 300, cfr = 0.5,
                          delay = Gamma(mean = 12.75, sd = 7))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  fit <- fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51),
                 cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0)
  expect_message(print(fit), "lognormal delay")
  expect_message(print(fit), "naive CFR")
  expect_invisible(suppressMessages(print(fit)))
})

test_that("a young outbreak is flagged low-information and print warns", {
  skip_if_no_cmdstan()
  set.seed(10)
  # A few very recent cases, almost nothing resolved: cfr stays near its prior.
  ll <- simulate_linelist(n = 25, cfr = 0.5, onset_days = 10,
                          delay = Gamma(mean = 12.75, sd = 7))
  d <- prepare_cfr_data(ll, obs_time = min(ll$onset_date) + 3)
  fit <- fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51),
                 cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0, seed = 1)
  expect_true(attr(summary(fit), "cfr_low_information"))
  expect_message(print(fit), "weakly identified")
})

test_that("a retrospective fit matches the naive proportion", {
  skip_if_no_cmdstan()
  set.seed(11)
  ll <- simulate_linelist(n = 600, cfr = 0.4,
                          delay = Gamma(mean = 12.75, sd = 7))
  d <- prepare_cfr_data(ll, obs_time = NULL)   # all resolved
  fit <- fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51),
                 cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0, seed = 1)
  s <- summary(fit)
  naive <- attr(s, "naive_cfr")
  # with every case resolved the cure model reduces to deaths / cases
  expect_equal(s[s$quantity == "cfr", "q50"], naive, tolerance = 0.02)
})

test_that("a real-time fit recovers the true CFR within its interval", {
  skip_if_no_cmdstan()
  set.seed(12)
  ll <- simulate_linelist(n = 800, cfr = 0.6, onset_days = 40,
                          delay = Gamma(mean = 12.75, sd = 7))
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  otd <- LogNormal(meanlog = Normal(2.41, 0.2),
                             sdlog = Normal(0.51, 0.15))
  fit <- fit_cfr(d, delay = otd, cfr_prior = Beta(1, 1),
                 chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0, seed = 1)
  cfr <- summary(fit)[1, ]
  expect_lt(cfr[["q2.5"]], 0.6)      # 95% CrI covers the true CFR
  expect_gt(cfr[["q97.5"]], 0.6)
})

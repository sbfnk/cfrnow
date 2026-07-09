skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  ok <- tryCatch(!is.null(cmdstanr::cmdstan_version()), error = function(e) FALSE)
  testthat::skip_if(!isTRUE(ok), "cmdstan not installed")
}

# shared sampler settings for the (slow) end-to-end fits
fit_quick <- function(d, ...) {
  fit_cfr(d, backend = "cmdstanr", chains = 2, iter = 800, warmup = 400,
          refresh = 0, seed = 1, ...)
}

test_that("a retrospective fit matches the naive proportion", {
  skip_if_no_cmdstan()
  set.seed(11)
  ll <- simulate_linelist(n = 1500, cfr = 0.4, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = NULL)          # every case resolved
  s <- summary(fit_quick(d))
  naive <- attr(s, "naive_cfr")
  # with every case resolved the cure model reduces to deaths / cases
  expect_equal(s$q50[s$quantity == "cfr"], naive, tolerance = 0.03)
})

test_that("real-time correction lifts the estimate above naive and covers truth", {
  skip_if_no_cmdstan()
  set.seed(12)
  ll <- simulate_linelist(n = 2000, cfr = 0.6, onset_days = 40,
                          delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  s <- summary(fit_quick(d))
  naive <- attr(s, "naive_cfr")
  cfr <- s[s$quantity == "cfr", ]
  expect_gt(cfr[["q50"]], naive)            # corrected > downward-biased naive
  expect_lt(cfr[["q2.5"]], 0.6)             # 95% CrI covers the true CFR
  expect_gt(cfr[["q97.5"]], 0.6)
  expect_true(all(c("rhat", "ess_bulk") %in% names(s)))
  expect_true(all(s$rhat < 1.05))
})

test_that("a gamma delay recovers the CFR and delay moments", {
  skip_if_no_cmdstan()
  set.seed(2)
  ll <- simulate_linelist(n = 2000, cfr = 0.4, onset_days = 40,
                          delay = Gamma(mean = 8, sd = 4))
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 2)
  s <- summary(fit_quick(d, family = stats::Gamma(link = "log")))
  expect_lt(s[s$quantity == "cfr", "q2.5"], 0.4)
  expect_gt(s[s$quantity == "cfr", "q97.5"], 0.4)
  expect_equal(s[s$quantity == "delay_mean", "q50"], 8, tolerance = 0.6)
})

test_that("the formula interface puts covariates on the CFR", {
  skip_if_no_cmdstan()
  set.seed(20)
  a <- simulate_linelist(n = 1200, cfr = 0.2, delay = LogNormal(2.4, 0.5))
  b <- simulate_linelist(n = 1200, cfr = 0.6, delay = LogNormal(2.4, 0.5))
  da <- as_epidist_cure_model(prepare_cfr_data(a, obs_time = NULL))
  db <- as_epidist_cure_model(prepare_cfr_data(b, obs_time = NULL))
  da$grp <- "low"; db$grp <- "high"
  d <- as_epidist_cure_model(rbind(da, db))
  fit <- fit_quick(d, formula = brms::bf(mu ~ 1, cfr ~ grp))
  fe <- brms::fixef(fit)
  expect_true("cfr_grplow" %in% rownames(fe))
  expect_lt(fe["cfr_grplow", "Estimate"], 0)          # low group < high group
})

test_that("a young outbreak is flagged low-information and print warns", {
  skip_if_no_cmdstan()
  set.seed(10)
  ll <- simulate_linelist(n = 25, cfr = 0.5, onset_days = 10,
                          delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = min(ll$onset_date) + 3)
  fit <- fit_quick(d)
  expect_true(attr(summary(fit), "cfr_low_information"))
  expect_message(print(fit), "weakly identified")
})

test_that("print reports the delay family, counts and naive CFR", {
  skip_if_no_cmdstan()
  set.seed(6)
  ll <- simulate_linelist(n = 400, cfr = 0.5, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = NULL)
  fit <- fit_quick(d)
  expect_message(print(fit), "lognormal delay")
  expect_message(print(fit), "naive CFR")
  expect_invisible(suppressMessages(print(fit)))
})

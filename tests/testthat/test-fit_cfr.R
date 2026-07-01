skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  ok <- tryCatch(!is.null(cmdstanr::cmdstan_version()), error = function(e) FALSE)
  testthat::skip_if(!isTRUE(ok), "cmdstan not installed")
}

test_that("retrospective fit recovers a known CFR", {
  skip_if_no_cmdstan()
  set.seed(1)
  ll <- simulate_linelist(n = 600, cfr = 0.45)
  d <- prepare_cfr_data(ll, obs_time = NULL)
  fit <- fit_cfr(d, chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0)
  s <- summarise_cfr(fit)
  cfr_row <- s[s$quantity == "cfr", ]
  expect_gt(cfr_row[["q97.5"]], 0.45)
  expect_lt(cfr_row[["q2.5"]], 0.45)
})

test_that("real-time correction lifts the estimate above the naive ratio", {
  skip_if_no_cmdstan()
  set.seed(2)
  # Growing epidemic sampled mid-outbreak: many recent, not-yet-resolved cases.
  ll <- simulate_linelist(n = 800, cfr = 0.6, onset_days = 40)
  cut <- max(ll$onset_date) - 2          # cut-off soon after the last onsets
  d <- prepare_cfr_data(ll, obs_time = cut)
  naive <- d$n_deaths / d$n_cases
  fit <- fit_cfr(d, chains = 2, parallel_chains = 2,
                 iter_warmup = 500, iter_sampling = 500, refresh = 0)
  s <- summarise_cfr(fit)
  cfr_med <- s[s$quantity == "cfr", "q50"]
  expect_gt(cfr_med, naive)             # corrected > downward-biased naive
})

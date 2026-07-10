test_that("pp_check_cfr rejects objects that are not cfrnow fits", {
  expect_error(pp_check_cfr(1), "must come from fit_cfr")
  expect_error(
    pp_check_cfr(structure(list(), class = "brmsfit")),
    "must come from fit_cfr"
  )
})

test_that("pp_check_cfr returns ggplots and reproduces the death count", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "cmdstan not installed"
  )

  set.seed(1)
  ll <- simulate_linelist(n = 300, cfr = 0.4, delay = LogNormal(2.4, 0.5))
  d <- prepare_cfr_data(ll, obs_time = max(ll$onset_date) - 5)
  fit <- fit_cfr(d,
    delay = LogNormal(meanlog = Normal(2.4, 0.2), sdlog = Normal(0.5, 0.15)),
    backend = "cmdstanr", chains = 1, iter = 400, warmup = 200,
    refresh = 0, seed = 1
  )

  # the cut-off is recorded so the check can replay the real-time truncation
  expect_false(is.na(fit$cfrnow$obs_time))
  expect_s3_class(pp_check_cfr(fit, "counts", ndraws = 50), "ggplot")
  expect_s3_class(pp_check_cfr(fit, "delay", ndraws = 50), "ggplot")

  # the observed death count should sit inside the posterior-predictive spread
  reps <- .cfr_replicate(fit, 200)
  pd <- reps$counts$n[reps$counts$outcome == "deaths"]
  obs <- reps$observed_counts$n[reps$observed_counts$outcome == "deaths"]
  expect_gte(obs, stats::quantile(pd, 0.01))
  expect_lte(obs, stats::quantile(pd, 0.99))
})

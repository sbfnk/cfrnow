# Input-validation guards and pure-R helpers. None of these need CmdStan, so
# they run everywhere and exercise the error paths taken before the sampler.

test_that("as_epidist_cure_model needs the model columns", {
  expect_error(as_epidist_cure_model(data.frame(x = 1)),
               "y|outcome|pwindow|swindow")
})

test_that("as_epidist_cure_model builds cure rows from prepare_cfr_data output", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    death_date = as.Date(c("2026-01-05", "2026-01-20", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  cure <- as_epidist_cure_model(d)
  expect_s3_class(cure, "epidist_cure_model")
  expect_equal(sum(cure$outcome == 1), d$n_deaths)   # 1 observed death
  expect_equal(sum(cure$outcome == 0), d$n_cens)     # future death + survivor
  expect_equal(nrow(cure), d$n_cases)
})

test_that("recovery-timing fits are refused on the epidist backend for now", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1),
    death_date = as.Date(c("2026-01-05", NA)),
    recovery_date = as.Date(c(NA, "2026-01-08"))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  expect_error(as_epidist_cure_model(d), "recovery-timing")
})

test_that("naive_cfr is deaths/cases, and NA when there are no cases", {
  expect_equal(naive_cfr(3, 10), 0.3)
  expect_true(is.na(naive_cfr(0, 0)))
})

test_that(".cfr_prior_sd reads a normal cfr prior and is NA otherwise", {
  expect_true(is.na(.cfr_prior_sd(brms::set_prior("beta(1, 1)", class = "b"))))
  s <- .cfr_prior_sd(cfr_default_prior())
  expect_true(is.finite(s) && s > 0 && s < 0.5)       # (0,1)-scale prior sd
})

test_that("simulate_linelist requires a delay", {
  expect_error(simulate_linelist(n = 5), "supply a `delay`")
})

test_that("prepare_cfr_data requires onset_date and death_date columns", {
  expect_error(prepare_cfr_data(data.frame(death_date = as.Date("2026-01-05"))),
               "onset_date")
  expect_error(prepare_cfr_data(data.frame(onset_date = as.Date("2026-01-01"))),
               "death_date")
})

test_that("summary.cfrnow_fit rejects objects not from fit_cfr", {
  expect_error(summary.cfrnow_fit(1), "fit_cfr")
})

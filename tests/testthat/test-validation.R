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

test_that("recovery dates produce timed-recovery rows flagged for two outcomes", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    death_date = as.Date(c("2026-01-05", NA, NA)),
    recovery_date = as.Date(c(NA, "2026-01-08", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  cure <- as_epidist_cure_model(d)
  expect_true(attr(cure, "use_recovery"))
  expect_equal(sum(cure$outcome == 2), d$n_recovery)   # timed recovery row
  expect_equal(sum(cure$outcome == 1), d$n_deaths)
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

test_that("cfr_prior moment-matches the CFR-scale mean and sd", {
  p <- cfr_prior(mean = 0.3, sd = 0.1)
  expect_equal(p$dpar, "cfr")
  expect_equal(.cfr_prior_sd(p), 0.1, tolerance = 0.02)   # induced sd, CFR scale
  m <- regmatches(p$prior,
                  regexec("normal\\(([^,]+), ([^)]+)\\)", p$prior))[[1]]
  set.seed(1)
  induced_mean <- mean(stats::plogis(stats::rnorm(1e5, as.numeric(m[2]),
                                                  as.numeric(m[3]))))
  expect_equal(induced_mean, 0.3, tolerance = 0.02)
})

test_that("fix_delay returns constant priors for the delay parameters", {
  p <- fix_delay(brms::lognormal(), mean = 12.75, sd = 7)
  expect_true(all(grepl("^constant\\(", p$prior)))
  expect_equal(sum(p$class == "Intercept"), 2L)
  expect_true(any(p$dpar == "sigma"))
  loc <- as.numeric(sub("constant\\(([^)]+)\\)", "\\1", p$prior[p$dpar == ""]))
  expect_equal(loc, log(12.75) - log1p((7 / 12.75)^2) / 2, tolerance = 1e-6)
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

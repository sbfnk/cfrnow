test_that("simulated delays match the requested mean and sd", {
  set.seed(1)
  specs <- list(LogNormal(mean = 12.75, sd = 7),
                Gamma(mean = 12.75, sd = 7))
  for (delay in specs) {
    ll <- simulate_linelist(n = 20000, cfr = 1, delay = delay)
    d <- as.numeric(ll$death_date - ll$onset_date)
    expect_equal(mean(d), 12.75, tolerance = 0.3)
    expect_equal(sd(d), 7, tolerance = 0.3)
  }
})

test_that("a recovery delay adds recoveries for non-fatal cases only", {
  set.seed(1)
  ll <- simulate_linelist(n = 300, cfr = 0.5,
                          delay = LogNormal(mean = 12.75, sd = 7),
                          recovery = LogNormal(mean = 21, sd = 9))
  expect_true("recovery_date" %in% names(ll))
  fatal <- !is.na(ll$death_date)
  expect_true(all(is.na(ll$recovery_date[fatal])))     # fatal: no recovery
  expect_true(all(!is.na(ll$recovery_date[!fatal])))   # non-fatal: recovered
  no_rec <- simulate_linelist(n = 5, delay = LogNormal(2, 0.5))
  expect_false("recovery_date" %in% names(no_rec))
})

test_that("unsupported families and prior-parameter delays are rejected", {
  expect_error(delay_native_order("weibull"))
  expect_error(stan_delay_fields(Normal(mean = 5, sd = 1),
                                 "dist_id", "p"))
  expect_error(sample_delay(5, LogNormal(
    meanlog = Normal(2, 0.1), sdlog = Normal(0.5, 0.1))))
})

test_that("stan_delay_fields reads estimated (Normal) native parameters", {
  d <- LogNormal(meanlog = Normal(2.4, 0.2),
                            sdlog = Normal(0.5, 0.15))
  ds <- stan_delay_fields(d, "dist_id", "p")
  expect_equal(ds$dist_id, primarycensored::pcd_stan_dist_id("lognormal", "delay"))
  expect_equal(ds$p1_est, 1L)
  expect_equal(ds$p2_est, 1L)
  expect_equal(ds$p1_prior_mean, 2.4)
  expect_equal(ds$p2_prior_sd, 0.15)
})

test_that("parse_cfr_prior accepts a fixed Beta and rejects everything else", {
  expect_equal(parse_cfr_prior(Beta(6.6, 13.4)),
               c(a = 6.6, b = 13.4))
  expect_equal(parse_cfr_prior(Beta(mean = 0.1, sd = 0.1))[["a"]],
               0.1 * (0.1 * 0.9 / 0.1^2 - 1))
  expect_error(parse_cfr_prior(5), "distspec distribution")
  expect_error(parse_cfr_prior(Gamma(shape = 3, rate = 1)), "Beta")
  expect_error(
    parse_cfr_prior(Beta(shape1 = Normal(2, 0.1),
                                    shape2 = 5)),
    "fixed"
  )
})

test_that("stan_delay_fields reads a fixed delay (fixed-F), with a q prefix too", {
  d <- Gamma(shape = 3, rate = 0.24)
  ds <- stan_delay_fields(d, "recovery_dist_id", "q")
  expect_equal(ds$recovery_dist_id,
               primarycensored::pcd_stan_dist_id("gamma", "delay"))
  expect_equal(ds$q1_est, 0L)
  expect_equal(ds$q2_est, 0L)
  expect_equal(ds$q1_fixed, 3)
  expect_equal(ds$q2_fixed, 0.24)
})
